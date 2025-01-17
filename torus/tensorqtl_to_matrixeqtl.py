import argparse
from collections import defaultdict
from fastparquet import ParquetFile
from glob import glob
import pandas as pd
# import subprocess

p = argparse.ArgumentParser(description="Convert gene-variant pairs from tensorQTL to Matrix eQTL output format.")
p.add_argument("permfile", help="tensorQTL top association per gene-variant pair file (*.cis_qtl.txt.gz). Used to subset data to eGenes.")
p.add_argument("parquets", help="Directory containing parquet files with all tested gene-variant pairs")
p.add_argument("outfile", help="Output file name (*.txt.gz)")
# p.add_argument("tmpfile1", help="Name of intermediate file to contain eQTL results in Matrix eQTL format")
# p.add_argument("tmpfile2", hep="Name of intermediate file to contain PIP values")
# p.add_argument("tissue", help="Name of tissue (will be included in output).")
args = p.parse_args()

# permfile = "~/br/data/tensorqtl/AQCT.cis_qtl.txt.gz"
# infile = "/home/dmunro/br/data/tensorqtl/AQCT/AQCT.cis_qtl_pairs.15.parquet"
# outfile = "~/br/data/dap/old/torus/Acbc.txt.gz"

perm = pd.read_csv(args.permfile, sep="\t")
perm = perm.loc[perm["qval"] < 0.05, :]
genes = set(perm["phenotype_id"])

df = []
for f in glob(f"{args.parquets}/*.parquet"):
    d = ParquetFile(f).to_pandas()
    d = d.loc[d["phenotype_id"].isin(genes)]
    # t-statistic is coefficient divided by its standard error (https://dss.princeton.edu/online_help/analysis/interpreting_regression.htm)
    d["t-stat"] = d["slope"] / d["slope_se"]
    d = d.rename(columns={"variant_id": "SNP", "phenotype_id": "gene", "slope": "beta", "pval_nominal": "p-value"})
    d = d[["SNP", "gene", "beta", "t-stat", "p-value"]]
    df.append(d)
df = pd.concat(df)
df.to_csv(args.outfile, sep="\t", index=False)

# subprocess.run(["torus/src/torus", "-d", args.tmpfile1, "-dump_pip", args.tmpfile2])
