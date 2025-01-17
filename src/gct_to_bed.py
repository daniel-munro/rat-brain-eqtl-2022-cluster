import argparse
import pandas as pd
import numpy as np
from gtfparse import read_gtf
import scipy.stats as stats

parser = argparse.ArgumentParser(
    description="Convert GCT gene expression to BED format and log-transform for TensorQTL."
)
parser.add_argument("infile", help="GCT file.")
parser.add_argument("annotations", help="GTF annotation file.")
parser.add_argument("outfile", help="Output file name ('*.bed').")
parser.add_argument(
    "--log", action="store_true", help="whether to log-transform values."
)
parser.add_argument(
    "--qnorm", action="store_true", help="whether to quantile normalize values."
)
args = parser.parse_args()

PSEUDOCOUNT = 1

# anno = read_gtf("data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf")
anno = read_gtf(args.annotations)
anno = anno[anno["feature"] == "gene"]
anno["start"][anno["strand"] == "-"] = anno["end"][anno["strand"] == "-"]
anno["end"] = anno["start"] + 1
anno = anno[anno["seqname"].isin([str(i) for i in range(1, 21)])]
anno["seqname"] = anno["seqname"].astype(int)  # for sorting
# anno = {
#     name: (chrom, tss)
#     for (name, chrom, tss)
#     in zip(anno["gene_name"], anno["seqname"], anno["start"])
# }
anno = anno.rename(columns={"seqname": "#chr"})
anno = anno[["#chr", "start", "end", "gene_id"]]
anno = anno.sort_values(["#chr", "start"])

# df = pd.read_csv("data/expression/ensemblGeneV2_Norm-counts-log2+1.txt", sep="\t")
df = pd.read_csv(args.infile, sep="\t", skiprows=2)
df.columns = df.columns.str.replace("_.*$", "")  # Remove region to match rat IDs in VCF
df = df.rename(columns={"Name": "gene_id"})
df = df.drop(columns="Description")

if args.qnorm:
    # From https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/rnaseqnorm.py
    R = stats.mstats.rankdata(df.iloc[:, 1:], axis=1)  # ties are averaged
    df.iloc[:, 1:] = stats.norm.ppf(R / (df.iloc[:, 1:].shape[1] + 1))
elif args.log:
    df.iloc[:, 1:] = np.log2(df.iloc[:, 1:] + PSEUDOCOUNT)

df = anno.merge(df, on="gene_id", how="inner")

# df.to_csv("data/expression/ensemblGeneV2_Norm-counts-log2+1.bed.gz", sep="\t",
#           index=False, float_format="%g")
df.to_csv(args.outfile, sep="\t", index=False, float_format="%g")
