import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Get top eQTL per gene from all chromosomes")
parser.add_argument("files", nargs="+", help="GEMMA output files")
parser.add_argument("-p", "--pvalue_col", required=False, help="Name of column with p-values used to choose top eQTLs. If omitted, keep all.")
parser.add_argument("-o", "--output", help="Output file (TSV)")
args = parser.parse_args()

top = []
for fname in args.files:
    # round_trip prevents precision issues in output
    d = pd.read_csv(fname, sep="\t", float_precision='round_trip')
    if args.pvalue_col is not None:
        d = d.groupby("gene_id")
        d = d.apply(lambda x: x.loc[x[args.pvalue_col] == min(x[args.pvalue_col]), :])
    # else:
        # d = d.loc[d[args.pvalue_col] < args.threshold, :]
    top.append(d)

pd.concat(top).to_csv(args.output, sep="\t", index=False)
