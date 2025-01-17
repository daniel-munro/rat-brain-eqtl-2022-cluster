import argparse
from fastparquet import ParquetFile

parser = argparse.ArgumentParser(
    description="Extract all pairs for a gene from full tensorQTL cis output."
)
parser.add_argument("parquet", help="Nominal parquet file containing gene.")
parser.add_argument("gene", help="ID of gene to extract.")
parser.add_argument("output", help="output file (TSV, can be *.gz)")
args = parser.parse_args()

sig = []
# Read each nominal file and save significant pairs.
nom = []
d = ParquetFile(args.parquet).to_pandas()
d = d.loc[d.phenotype_id == args.gene, :]

d.to_csv(args.output, sep="\t", index=False, float_format="%g")
