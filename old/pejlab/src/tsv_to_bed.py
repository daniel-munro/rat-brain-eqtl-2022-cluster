import argparse
import pandas as pd
import numpy as np
from gtfparse import read_gtf

parser = argparse.ArgumentParser(description="Convert TSV gene expression to BED format for TensorQTL.")
parser.add_argument("infile", help="TSV file. Row names are gene IDs, so there is 1 less column ID than columns.")
parser.add_argument("annotations", help="GTF annotation file.")
parser.add_argument("outfile", help="Output file name ('*.bed').")
parser.add_argument("--region", help="Name of brain region to subset to.")
args = parser.parse_args()

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
df = pd.read_csv(args.infile, sep="\t")
df = df.filter(like=args.region, axis=1)
df.columns = df.columns.str.replace("_.*$", "")
# Currently TensorQTL fails when there are constant phenotypes AND other genes with empty cis-windows:
# df = df[np.sum(df != 0, axis=1) > 4]

df.index.name = "gene_id"
df = anno.merge(df, on="gene_id", how="inner")

# df.to_csv("data/expression/ensemblGeneV2_Norm-counts-log2+1.bed.gz", sep="\t",
#           index=False, float_format="%g")
df.to_csv(args.outfile, sep="\t", index=False, float_format="%g")
