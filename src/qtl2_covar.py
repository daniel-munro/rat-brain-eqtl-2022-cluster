import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]
columns = sys.argv[3].split(",")

df = pd.read_csv(infile)
df = df.rename(columns={"rat_rfid": "id", "batchnumber": "batch", "sex_mf": "sex"})
df["id"] = df["id"].str.upper()
# df = df[["id", "batch", "sex"]]
df = df[["id"] + columns]
df["generations"] = 90
df.to_csv(outfile, index=False)
