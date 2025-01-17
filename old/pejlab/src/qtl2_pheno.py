import sys
import pandas as pd
import numpy as np

infile = sys.argv[1]
outfile = sys.argv[2]
pseudocount = 5

df = pd.read_table(infile, skiprows=2)
df = df.rename(columns={"Name": "id"})
df = df.drop(columns="Description")
for c in df.columns[1:]:
    # if df[c].dtype == "int64":
    df[c] = np.log2(pseudocount + df[c])
    df = df.rename(columns={c: c.split("_")[0]}) # Must be rat ID to match geno file.
df.to_csv(outfile, index=False, float_format="%.6f")
