import os
import sys
import argparse
import pandas as pd
import numpy as np
import scipy
import scipy.interpolate
import pickle

# Down sample the data by averaging the value in size r window
def down_sample(x, r):
    # x: input 2D numpy array
    # r: ratio of down sampling
    ns = x.shape[1] // r       # number of samples
    #print(f"ns = {ns}")
    # Remove remainder columns more than ns * r
    x = x[:, 0:int(ns*r)]
    # Average every r elements
    x = x.reshape(x.shape[0], -1, r).mean(axis=2)
    return x


##
def read_gene_structs(bed_file="Ehux_genbank.bed"):
    # Read the gene structures from the input bed file
    genes = pd.read_csv(bed_file, sep="\t", header=None, index_col=3)
    genes.columns = ["chr", "start", "end", "name", "strand", "start_g", "end_g", "rgb", "num_exons", "exon_length",
                     "exon_start"]
    # Read the mapping from genebank rna names to E. hux gene IDs
    mapping = pd.read_csv("genbank_mapping.txt", header=None, sep="\t", index_col=0)

    # some gene IDs are duplicated because of alternative splicing.
    # In this case, we only keep the first one
    seen_id = set()
    uniq_rna = []
    dup_rna = []
    for x in mapping.index:
        if x in genes.index:
            if mapping.loc[x, 2] not in seen_id:
                uniq_rna.append(x)
                seen_id.add(mapping.loc[x, 2])
            else:
                dup_rna.append(x)
    genes = genes.loc[uniq_rna]

    #print(mapping.loc[dup_rna])
    # mapping rna names to Ehux gene IDs
    genes["id"] = mapping.loc[genes.index, 2]
    genes["rna"] = genes.index
    genes.index = genes["id"]
    # Extract only the useful columns
    genes = genes[["rna", "chr", "start", "end", "strand", "num_exons", "exon_length", "exon_start"]]

    return genes


# Return the position of Transcription starting site (TSS)
# The TSS position depends on whether the gene is on the pos or neg strand
def tss_pos(gene):
    if gene.strand == '+':
        return gene.start
    else:
        return gene.end


# Calculate the start and the end of the window flanking the TSSs
# Only choose the genes whose windows are within the chr boundary
def calc_tss_windows(genes, win_size):
    window_start = [tss_pos(row) - win_size for index, row in genes.iterrows()]
    window_end = [tss_pos(row) + win_size for index, row in genes.iterrows()]
    genes["window_start"] = window_start
    genes["window_end"] = window_end
    return genes


# Filter genes based on selected gene IDs from a file
def filter_genes_by_ids(genes, sel_id_file):
    # Read the file of selected gene IDs
    sel_genes = pd.read_csv(sel_id_file, sep="\t", header=0, index_col=0)
    sel_ids = list(set(genes.index) & set(sel_genes.index))
    genes_sel = genes.loc[sel_ids]
    # Add a column of "log2FC"
    genes_sel["LFC"] = sel_genes.loc[sel_ids, "log2FC"]
    return genes_sel


# Filter genes whose windows are larger than the sizes of the chrs
def filter_genes_by_window(genes, win_size):
    # calculate window surrounding the gene TSSes
    calc_tss_windows(genes, win_size)
    # Read chr lengths from a file
    chr_lens = pd.read_csv("Ehux_genome_length.txt", sep="\t", header=None, index_col=0)
    # Filter out genes with window_start < 0 and window_end > corresponding chr length
    genes['chr_len'] = [chr_lens.loc[row.chr, 1] for index, row in genes.iterrows()]
    genes = genes[genes['window_start'] >= 0]
    genes = genes[genes['window_end'] < genes['chr_len']]
    return genes


def filter_genes_by_length(genes, min_length=1000):
    genes = genes[(genes["end"] - genes["start"]) >= min_length]
    return genes


# Read methylation data files 
# the methylation data is organized into dictionary with an entry for each chr
def read_methylation_data(methy_file):
    meth = {}
    # Read the file into a dataframe
    data = pd.read_csv(methy_file, header=0, sep="\t", index_col=None)

    while data.shape[0] > 0:
        cur_chr = data.iloc[0]["V1"]
        i = 0
        for index, row in data.iterrows():
            if row["V1"] != cur_chr:
                break
            i = i + 1
        print(cur_chr, i)
        meth[cur_chr] = data.iloc[0:i]
        # remove data that has been already processed
        data = data.iloc[i:]

    #chrs = set(data["V1"])
    #for chr in chrs:
    #    meth[chr] = data[data["V1"]== chr]
    """"
    for index, row in data.iterrows():
        if row.V1 not in meth:
            meth[row.V1] = pd.DataFrame(columns=data.columns).append(row)
        else:
            meth[row.V1] = meth[row.V1].append(row)
    """
    # Use the start position within chr as index
    for chr, data in meth.items():
        data.index = data["V2"]

    return meth

# Get interpolated methylation values in a window
# fm - dataframe of methylation data
# flank_win - extension of window in both directions for interpolation
# interp_method - method of interpolation, can be either "linear", "pchip" or "none"
# Return - interpolated methylation values between win_start and win_end
def interpolate_methyl_window(
        fm, win_start, win_end, flank_win=5000,
        interp_method="linear",
        slide_win_size=50,
):
    interp_win_start = max(0, win_start - flank_win)
    interp_win_end = win_end + flank_win
    raw_fm = fm.loc[interp_win_start:interp_win_end]
    raw_x = np.array(raw_fm.index)
    raw_y = np.array(raw_fm)

    interp_y = np.zeros(win_end - win_start)
    # Linear 1D interpolation
    if interp_method == "linear":
        interp_f = scipy.interpolate.interp1d(raw_x, raw_y)

        # the [win_start, win_end] window may be larger than the available data in raw_x and raw_y
        x = range(max(win_start, raw_x[0]), min(win_end, raw_x[-1]))
        y = interp_f(x)
        #print(f"x[0] = {x[0]}, x[-1] = {x[-1]}")
        interp_y[x[0] - win_start : x[-1] + 1 - win_start] = y
    elif interp_method == "pchip":
        x = range(win_start, win_end)
        interp_y = scipy.interpolate.pchip_interpolate(raw_x, raw_y, x)
    elif interp_method == "avg":  # calculate average methylation in sliding windows
        """fm_window = fm.loc[win_start:win_end]
        count_y = np.zeros(win_end - win_start)
        for pos, meth in fm_window.iteritems():
            idx = pos - win_start
            if idx < win_end - win_start:
                interp_y[idx] = meth
                count_y[idx] = 1
        interp_y = np.reshape(interp_y, (-1, slide_win_size)).sum(axis=1)
        count_y = np.reshape(count_y, (-1, slide_win_size)).sum(axis=1)
        count_y = np.clip(count_y, 1, None)  # set min count to 1 to avoid dividing by zero
        interp_y = interp_y / count_y
        interp_y = np.repeat(interp_y, slide_win_size)"""
        """for x in range(win_start, win_end):
            sw_fm = raw_fm.loc[x-slide_win_size:x+slide_win_size]
            count = max(1, sw_fm.shape[0])
            interp_y[x-win_start] = np.array(sw_fm).sum() / count"""
        raw_y = np.zeros(win_end - win_start + 2 * slide_win_size)
        cnt_y = np.zeros(win_end - win_start + 2 * slide_win_size)
        sw_fm = raw_fm.loc[(win_start - slide_win_size):(win_end+slide_win_size)]
        for pos, meth in sw_fm.iteritems():
            idx = pos - (win_start-slide_win_size)
            if idx < (win_end-win_start+2*slide_win_size):
                raw_y[idx] = meth
                cnt_y[idx] = 1
        for i in range(win_end-win_start):
            interp_y[i] = (raw_y[i:(i+2*slide_win_size)]).sum()
            count = max(1, (cnt_y[i:(i+2*slide_win_size)]).sum())
            interp_y[i] = interp_y[i]/count
    else: # no interpolation
        fm_window = fm.loc[win_start:win_end]
        for pos, meth in fm_window.iteritems():
            idx = pos - win_start
            interp_y[idx] = meth

    return interp_y


# Calculate methylation changes in a window surrounding TSS
def calc_methylation_diff_in_window(genes, min_cpgs=30, interp_method="linear",
                                    meth_in_file="Ehux_meth_merged.tsv",
                                    md_03_sites=3,):
    # Read methylation data files
    meth_pickle_file = meth_in_file+".pickle"
    if not os.path.exists(meth_pickle_file):  # Pickle file doesn't exist
        print(f"{meth_pickle_file} does not exist")
        #meth_1516 = read_methylation_data("EH1516_meth_merged.tsv")
        #meth_217 = read_methylation_data("EH217_meth_merged.tsv")
        meth = read_methylation_data(meth_in_file)
        # Save the meth_data because it takes a long time to compute
        #save_meth = {"meth_1516": meth_1516, "meth_217": meth_217}
        save_meth = {"meth": meth}
        pickle.dump(save_meth, open(meth_pickle_file, "wb"), protocol=3)
    else:
        with open(meth_pickle_file, "rb") as fp:
            save_meth = pickle.load(fp)
            #meth_1516 = save_meth["meth_1516"]
            #meth_217 = save_meth["meth_217"]
            meth = save_meth["meth"]
    # The total window size surrounding the TSS (window_end is not inclusive)
    win_len = genes.iloc[0]['window_end'] - genes.iloc[0]['window_start']
    num_gene = genes.shape[0]
    print(f"# of genes = {num_gene}, window size = {win_len}")
    #
    gene_meth_diff = {}
    for id, gene in genes.iterrows():
        vals = np.zeros(win_len)
        
        # Select methyl data within the window range
        if not (gene.chr in meth):     # if methylation data isn't available for the given chr, just continue
            continue
        fm = meth[gene.chr]

        # Select CpGs within the window range
        fm_window = fm.loc[gene.window_start : gene.window_end]
        # Filter out a window if the number of CpGs < min_meth
        print(f"gene {id} has {fm_window.shape[0]} CpG sites")
        #print(fm_window.index)
        if fm_window.shape[0] < min_cpgs:
            sys.stderr.write(
                f"CpG sites Warning: {id} has {fm_window.shape[0]} < {min_cpgs} CpGs\n"
            )
            continue

        # Filter out a window if there are fewer than 3 CpGs in the window with |meth_diff| > 0.3
        window_meth_diff = fm_window["beta.217"] - fm_window["beta.1516"]
        md_sites = (window_meth_diff.abs() >= 0.3).sum()
        if md_sites < md_03_sites:
            sys.stderr.write(
                f"CpG sites Warning: {id} has {md_sites} < {md_03_sites} CpGs with methylation changes > 0.3 \n"
            )
            continue

        #print(f"gene={id}, {gene.chr}, win_start={gene.window_start}, win_end={gene.window_end}")
        fm_win_1516 = interpolate_methyl_window(fm["beta.1516"], gene.window_start, gene.window_end, interp_method=interp_method)
        fm_win_217 = interpolate_methyl_window(fm["beta.217"], gene.window_start, gene.window_end, interp_method=interp_method)

        vals = fm_win_217 - fm_win_1516
        # Reverse the vals array for genes on the negative strand
        if gene.strand == '-':
            vals = vals[::-1]
        gene_meth_diff[id] = vals
    # Return results as a dataframe
    return pd.DataFrame.from_dict(gene_meth_diff, orient="index")

# randomly shuffle and divide the data into training, validation and test sets
def data_prep(meth_diff, expr_lfc):
    # Use a seed for reproducibility
    np.random.seed(42)
    # shuffle the dataset randomly
    permut = np.random.permutation(meth_diff.shape[0])
    print("permutations:", permut[0:10])
    meth_diff = meth_diff.iloc[permut]
    expr_lfc = expr_lfc.iloc[permut]
    print(expr_lfc.head())
    print(np.sign(expr_lfc.head()))
    print(meth_diff.head())

    # split data into training(60%), validation(20%) and testing(20%) sets
    l = meth_diff.shape[0]
    lv = int(l / 5)  # number of validation data
    lt = int(l / 5)  # number of testing data
    ln = l - lv - lt  # number of training data

    start, end = 0, ln
    print(f"start: {start}, end: {end}")
    mdf_train = np.asarray(meth_diff.iloc[start:end])
    lfc_train = np.asarray(expr_lfc.iloc[start:end])
    print(mdf_train.shape)
    print(mdf_train[0:5, 0:5])
    start, end = ln, ln + lv
    print(f"start: {start}, end: {end}")
    mdf_valid = np.asarray(meth_diff.iloc[start:end])
    lfc_valid = np.asarray(expr_lfc.iloc[start:end])
    start, end = ln + lv, ln + lv + lt
    print(f"start: {start}, end: {end}")
    mdf_test = np.asarray(meth_diff.iloc[start:end])
    lfc_test = np.asarray(expr_lfc.iloc[start:end])

    print(mdf_train.shape, lfc_train.shape)
    print(mdf_train[0:5, 0:5])
    print(np.min(mdf_train), np.max(mdf_train))
    print(np.min(lfc_train), np.max(lfc_train))

    return [mdf_train, lfc_train, mdf_valid, lfc_valid, mdf_test, lfc_test]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract methylation difference data from two E. hux samples"
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="in_file",
        default="Ehux_meth_merged.tsv",
        help="Input methylation difference data",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        default=None,
        help="Output file of methylation data and logFC divided into training, validation and testing datasets",
    )
    parser.add_argument(
        "--ids",
        dest="geneid_file",
        default=None,
        help="List of selected gene IDs",
    )
    parser.add_argument(
        "-w",
        "--window-size",
        dest="window_size",
        type=int,
        default=3000,
        help="Upstream and downstream window size around TSS",
    )
    parser.add_argument(
        "--down-sampling",
        dest="down_sampling",
        type=int,
        default=20,
        help="Ratio of down sampling the methylation data",
    )
    parser.add_argument(
        "--min-cpg-sites",
        dest="min_cpg_sites",
        type=int,
        default=30,
        help="Minimum number of cpg sites in the [-window_size, window_size] promoter region for genes to be considered",
    )
    parser.add_argument(
        "--interp-method",
        dest="interp_method",
        default="linear",
        help="Interpolation method of methylation data",
    )
    parser.add_argument(
        "--threshold_meth_diff",
        dest="threshold_meth_diff",
        type=float,
        default=0.3,
        help="",
    )

    args = parser.parse_args()

    # Read gene structures from a bed file
    genes = read_gene_structs()

    # Filter genes based on selected gene IDs from a file
    genes = filter_genes_by_ids(genes, args.geneid_file)

    # Filter genes based on the windows surrounding the TSSs
    genes = filter_genes_by_window(genes, args.window_size)
    print(genes.shape)

    # Calcualte methylation difference in gene windows
    meth_diff = calc_methylation_diff_in_window(genes, min_cpgs=args.min_cpg_sites,
                                                interp_method=args.interp_method,
                                                meth_in_file=args.in_file,
                                                md_03_sites=3,
                                                )
    meth_diff = meth_diff.sort_index()
    print(meth_diff.shape)

    # Filter out all rows that don't have significant methylation changes,
    # i.e. the max absolute methyl diff < 0.3
    #diff_methyled = meth_diff.abs().max(axis=1) >= args.threshold_meth_diff
    #meth_diff = meth_diff.loc[diff_methyled, :]
    #print(meth_diff.shape)

    # Only keep the genes that have methylation diff data
    expr_lfc = genes.loc[meth_diff.index, "LFC"]
    expr_lfc = expr_lfc.sort_index()
    print(expr_lfc.shape)
    print(expr_lfc.head(5))

    # Down sampling the methylation difference data
    if args.down_sampling > 1:
        data_ds = down_sample(meth_diff.to_numpy(), args.down_sampling)
        meth_diff = pd.DataFrame(data_ds, index=meth_diff.index)
        print(f"After down sampling, {meth_diff.shape}")
        print(meth_diff.iloc[0:5, 100:105])

    [mdf_train, lfc_train, mdf_valid, lfc_valid, mdf_test, lfc_test] = data_prep(meth_diff, expr_lfc)

    # Save the methylation diff and gene log2FC data into a pickle file
    if(args.out_file == None):
        args.out_file = "data_EH1516_EH217_"+args.interp_method+".pickle"
    #pickle_file = "data_EH1516_EH217_filtered.pickle"
    with open(args.out_file, "wb") as fp:
        #save = {"meth_diff": meth_diff, "expr_lfc": expr_lfc}
        save = {
            "meth_diff": meth_diff,
            "expr_lfc": expr_lfc,
            "train_dataset": mdf_train,
            "train_labels": lfc_train,
            "valid_dataset": mdf_valid,
            "valid_labels": lfc_valid,
            "test_dataset": mdf_test,
            "test_labels": lfc_test,
        }
        pickle.dump(save, fp)
        # pickle.dump(save, f)
        fp.close()




