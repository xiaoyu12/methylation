from extract_meth_diff_data import read_gene_structs
from extract_meth_diff_data import filter_genes_by_ids, filter_genes_by_length, filter_genes_by_window

def is_overlap(gene1, gene2):
    if gene1.chr != gene2.chr or gene1.strand != gene2.strand:
        return False
    if gene1.start > gene2.end:
        return False
    if gene2.start > gene1.end:
        return False
    return True


if __name__ == '__main__':
    # Read gene structures from a bed file
    genes = read_gene_structs()

    # Filter genes based on selected gene IDs from a file
    genes = filter_genes_by_ids(genes, "DE_E1516_vs_E217.annot.txt")

    print(genes.shape)

    # Filter genes based on the windows surrounding the TSSs
    genes = filter_genes_by_window(genes, 3000)
    print(genes.shape)

    # Filter out genes with length < 1000
    genes = filter_genes_by_length(genes, min_length=1000)
    print(genes.shape)

    # Add genes to a dictionary
    gene_dict = {}
    for id, g in genes.iterrows():
        #print(g.chr)
        if g.chr not in gene_dict:
            gene_dict[g.chr] = [g]
        else:
            gene_dict[g.chr].append(g)
            #print(len(gene_dict[g.chr]))

    # Find if there are any overlaps among the genes
    overlapped = []
    overlapped_ids = []
    for key, list in gene_dict.items():
        #print(key)
        for i in range(len(list)):
            for j in range(i+1, len(list)):
                if is_overlap(list[i], list[j]):
                    #print("overlap", list[i], list[j])
                    overlapped.append(list[i])
                    overlapped.append(list[j])

    for g in overlapped:
        print(g.name)
        overlapped_ids.append(g.name)
    
    # Remove genes that have overlaps 
    sel_ids = set(genes.index).difference(set(overlapped_ids))
    genes = genes.loc[sel_ids]
    print(genes.shape)
    
    genes.to_csv("DE_no_overlap.tsv", sep="\t")
    """for i in range(genes.shape[0]):
        gene1 = genes.iloc[i]
        print(i)
        for j in range(i+1, genes.shape[0]):
            gene2 = genes.iloc[j]
            if is_overlap(genes.iloc[i], genes.iloc[j]):
                print ("overlap", gene1, gene2)
                """