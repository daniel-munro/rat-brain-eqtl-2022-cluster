import argparse
import pandas as pd
import numpy as np
from gtfparse import read_gtf

parser = argparse.ArgumentParser(
    description="Convert TSV gene expression to BED format for TensorQTL."
)
parser.add_argument(
    "infile",
    help="TSV file. Row names are gene IDs, so there is 1 less column ID than columns.",
)
parser.add_argument("annotations", help="GTF annotation file.")
parser.add_argument("outfile", help="Output file name ('*.bed').")
args = parser.parse_args()

anno = read_gtf(args.annotations)
anno = anno[anno["feature"] == "gene"]
anno["start"][anno["strand"] == "-"] = anno["end"][anno["strand"] == "-"]
anno["end"] = anno["start"] + 1
anno = anno[anno["seqname"].isin([str(i) for i in range(1, 21)])]
anno["seqname"] = anno["seqname"].astype(int)  # for sorting
anno = anno.rename(columns={"seqname": "#chr"})
anno = anno[["#chr", "start", "end", "gene_id"]]
anno = anno.sort_values(["#chr", "start"])

df = pd.read_csv(args.infile, sep="\t")

## This has already been done for the original expression files.
# ## Address sample mixups with swap/removal
# if np.all(np.isin(["000789FFF0_LHB", "000789FFF9_LHB"], df.columns)):
#     tmp = df["000789FFF0_LHB"].copy()
#     df["000789FFF0_LHB"] = df["000789FFF9_LHB"]
#     df["000789FFF9_LHB"] = tmp
# remove = ["00078A0224_Acbc", "000789FFF8_IL", "00078A2667_Acbc", "00078A2667_IL",
#           "00078A2667_PL", "00078A2667_VoLo", "00078A18A7_Acbc", "00078A18A7_IL",
#           "00078A18A7_LHB", "00078A18A7_PL", "00078A18A7_VoLo"]
# df = df.drop(remove, axis=1, errors="ignore")

df.columns = df.columns.str.replace("_.*$", "", regex=True)
# Currently TensorQTL fails when there are constant phenotypes AND other genes with empty cis-windows:
# TensorQTL perm beta parameter fitting failed on a gene with only 1 nonzero value.
df = df.loc[np.sum(df != 0, axis=1) > 1, ]

df.index.name = "gene_id"
df = anno.merge(df, on="gene_id", how="inner")

df.to_csv(args.outfile, sep="\t", index=False, float_format="%g")
