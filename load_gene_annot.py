import re
import sys
import os
import argparse
import pandas as pd
import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF

# import gffutils
import gzip


class Transcript:
    def __init__(
        self,
        transcript_id,
        chrom,
        strand,
        start,
        end,
        exons,
        cds,
        cds_tx_start_cmpl,
        cds_tx_end_cmpl,
        gene_id=None,
        common_name="",
    ):
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end
        self.exons = exons
        self.cds = cds
        self.cds_tx_start_cmpl = cds_tx_start_cmpl
        self.cds_tx_end_cmpl = cds_tx_end_cmpl
        self.gene_id = gene_id
        self.common_name = common_name

    # Str for print a transcript
    def __str__(self):
        ex_starts = ",".join([str(x.start) for x in self.exons])
        ex_ends = ",".join([str(x.end) for x in self.exons])
        s = "\t".join(
            [
                str(x)
                for x in [
                    self.transcript_id,
                    self.gene_id,
                    self.chrom,
                    self.start,
                    self.end,
                    self.strand,
                    self.end - self.start + 1,
                    len(self.exons),
                    ex_starts,
                    ex_ends,
                ]
            ]
        )
        return s


class Gene:
    def __init__(self, gene_id):
        self.gene_id = gene_id
        # A gene may have multiple transcripts
        self.transcripts = list()

    # Add a transcript to the list
    def append_transcript(self, transcript):
        self.transcripts.append(transcript)

    #
    def __str__(self):
        s = ""
        for transcript in self.transcripts:
            s = s + str(transcript)
        return s


class TranscriptFilter:
    def __init__(
        self,
        only_simple_chrs=True,
        only_mrna=True,
        only_complete=True,
        only_uniq_tx_start=False,
        collapse_on_common_name=False,
        exclude_lt_four_exons=False,
        gene_list_only=None,
    ):
        ### Transcript Level Filters
        ## Only look at transcripts defined as complete. cds start and end well defined.
        self.only_complete = only_complete
        ### Gene Level Filters
        ## Only look at genes on chromosomes defined in valid_chr
        self.only_simple_chrs = only_simple_chrs
        ## Only look at genes defined as protein coding (starting with NM_)
        self.only_mrna = only_mrna
        ## Only look at genes with unique TxStart sites
        self.only_uniq_tx_start = only_uniq_tx_start
        ## Collapse Annotations so that only common gene name represented
        self.collapse_on_common_name = collapse_on_common_name
        ## Exclude genes with < 4 exons to have straightforward comparison with ROI
        self.exclude_lt_four_exons = exclude_lt_four_exons
        ## Only consider genes in a given list
        self.gene_list_only = gene_list_only

    # Check if a gene_id is in the given list
    def filter_gene_list(self, gene_id):
        if self.gene_list_only is None:  # No filtering for gene list
            return True
        return gene_id in self.gene_list_only

    # Check whether a transcript passes the filtering
    def filter_transcript(self, transcript):
        if self.only_mrna:  # Check if the transcript is protein coding
            if not self.transcript_is_coding(transcript):
                return False
        if (
            self.only_complete
        ):  # Check if the transcript has annotated start and stop codon
            if not self.transcript_is_complete(transcript):
                return False
        if self.only_simple_chrs:
            if not self.transcript_on_simple_chrom(transcript):
                return False
        return True

    # Check if a transcript is protein coding
    def transcript_is_coding(self, transcript):
        return len(transcript.cds) > 0

    # Check if a transcript has annotated start and stop codon
    def transcript_is_complete(self, transcript):
        return transcript.cds_tx_start_cmpl and transcript.cds_tx_end_cmpl

    # Check if a transcript is located on a simple chromosome
    def transcript_on_simple_chrom(self, transcript):
        return transcript.chrom in valid_chr


class GeneAnnotLoader:
    def __init__(self, transcript_filter=TranscriptFilter()):
        # Default flags to filter gene annotations
        self.transcript_filter = transcript_filter
        self.genes = {}

    # Read gene annotations from the input file
    def load(self, annot_file):
        sys.stderr.write("Processing Annotation from %s\n" % (annot_file))
        # Count the number of lines in the annotation file for each chromosome
        chr_set = {}
        with open(annot_file) as fh:
            for line in fh:
                chr = line.strip().split()[0]
                if not (chr in chr_set.keys()):
                    chr_set[chr] = 1
                else:
                    chr_set[chr] += 1
        print(len(chr_set))
        for key, count in chr_set.items():
            print(key, count)
        # Parse the input GTF/GFF file
        self.parse_gtf_file(annot_file)
        print("# of parsed genes = ", len(self.genes))
        # for k, v in self.genes.items():
        #    print("gene ", k)
        #    print(v)
        return self.genes

    ## Parse the input GTF file
    def parse_gtf_file(self, annot_file):
        # examiner = GFFExaminer()
        # fh = open(annot_file)
        # pprint.pprint(examiner.available_limits(fh))
        # fh.close()
        count = 0
        with open(annot_file) as fh:
            # for rec in GFF.parse(fh, target_lines=1000):
            for rec in GFF.parse(fh):
                self.chrom = rec.id
                # for rec in GFF.parse(fh):
                print("Processing chrom ", rec.id)
                for feature in rec.features:
                    try:
                        transcript = self.transcript_from_gtf_feature(feature)
                        # print(transcript)
                        gene_id = transcript.gene_id
                        if self.transcript_filter.filter_gene_list(gene_id):
                            # print(gene_id, "is in the list")
                            # Add the gene if it isn't in the list yet
                            if gene_id not in self.genes.keys():
                                self.genes[gene_id] = Gene(gene_id)
                            # Append the transcript to its parent gene
                            self.genes[gene_id].append_transcript(transcript)
                            # print("# genes parsed: ", len(self.genes))
                        else:
                            # print(gene_id, "not found in the list")
                            continue
                    except:
                        continue

        # db = gffutils.create_db(annot_file, ":memory:")

    ## Construct a Transcript object from a parsed GTF feature
    def transcript_from_gtf_feature(self, feature):
        # transcript id
        tx_id = feature.id
        if feature.sub_features[0].type != "transcript":  # not a transcript feaure
            return None
        gene_id = feature.sub_features[0].qualifiers["gene_id"][0]
        strand = feature.strand
        exons = []
        cds = []
        pos_low = feature.location.start
        pos_high = feature.location.end
        cds_tx_start_cmpl = False
        cds_tx_end_cmpl = False
        for feat in feature.sub_features:
            if feat.type == "exon":
                exons.append(feat.location)
            elif feat.type == "CDS":
                cds.append(feat.location)
            elif feat.type == "start_codon":
                cds_tx_start_cmpl = True
            elif feat.type == "stop_codon":
                cds_tx_end_cmpl = True

        transcript = Transcript(
            tx_id,
            self.chrom,
            strand,
            pos_low,
            pos_high,
            exons,
            cds,
            cds_tx_start_cmpl,
            cds_tx_end_cmpl,
            gene_id=gene_id,
        )
        return transcript

    ## Filtering the parsed gene list
    def filter_genes(self):
        valid_genes = {}
        for gene_id, gene in self.genes.items():
            for transcript in gene.transcripts:
                if self.transcript_filter.filter_transcript(
                    transcript
                ):  # transcript is valid
                    valid_genes[
                        gene_id
                    ] = transcript  # the first transcript passed the filtering
                    break
        return valid_genes

    ## Read lines from a GTF file and output only lines for gene_ids in a given list
    ## The purpose of this function to reduce the GTF file size for further analysis
    def filter_gtf_only_in_genlist(self, in_gtf_file, out_gtf_file):
        with open(out_gtf_file, "w") as output:
            with open(in_gtf_file) as input:
                for line in input:
                    gene_id = line.strip().split()[9]
                    # m = re.match(r"gene_id \"(\w+)\";")
                    gene_id = gene_id[1:-2]
                    if self.transcript_filter.filter_gene_list(gene_id):
                        output.write(line)


def setup_parser_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Load methylation and gene expression data",
    )
    # parser.add_argument("bedgraph", help="Sample bedgraph or DMR (DSS) output file")
    parser.add_argument(
        "-g",
        "--gene_annot",
        dest="annot_file",
        help="Gene annoation file (REFSEQ or ENSG GTF format)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        default="./out.txt",
        help="Output file of selected genes and transcripts",
    )
    return parser


### Define valid chromosomes
valid_chr = [
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
    ]

# main function
if __name__ == "__main__":
    parser = setup_parser_arguments()
    args = parser.parse_args()

    # Read the expression file of protein coding genes and only consider genes in this list
    rpkm = pd.read_table(
        "data/expression/57epigenomes.RPKM.pc.tsv", sep="\t", header=0, index_col=0
    )
    gene_ids = set(rpkm.index)
    print("num of genes: " + str(len(gene_ids)))
    filter_flags = TranscriptFilter(gene_list_only=gene_ids)
    # filter_flags = FilterFlags()

    # Load gene annotations from a file (GTF)
    gene_annot_loader = GeneAnnotLoader(transcript_filter=filter_flags)
    gene_annot_loader.load(args.annot_file)
    valid_genes = gene_annot_loader.filter_genes()
    print("# of genes after filtering: ", len(valid_genes))
    # Print out select genes and transcripts to a file
    with open(args.out_file, "w") as outfile:
        header_list = [
            "TRANSCRIPT_ID",
            "GENE_ID",
            "CHR",
            "START",
            "END",
            "STRAND",
            "LENGTH",
            "NUM_EXONS",
            "EXON_STARTS",
            "EXON_ENDS",
        ]
        outfile.write("#" + "\t".join(header_list) + "\n")
        for gid, transcript in valid_genes.items():
            outfile.write(str(transcript) + "\n")

    # gene_annot_loader.filter_gtf_only_in_genlist(args.annot_file, "data/hg19.ensGene.pc.gtf")
