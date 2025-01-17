import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

# Was using DemographicInfo_TN_RNASeq_88_LauraSaba.csv, but switched to rat_info.txt.
df = pd.read_csv(infile, sep="\t")

# df = df.rename(columns={"rat_rfid": "ID", "batchnumber": "batch", "sex_mf": "sex"})
# df["ID"] = df["ID"].str.upper()  # Some letters lowercase for some reason.
# df = df[["ID", "batch", "sex"]]
# df = pd.get_dummies(df, columns=["batch", "sex"], drop_first=True)

df = df.rename(columns={"rat_id": "ID", "sex_mf": "sex"})
df = df[["ID", "sequencing_batch", "rat_batch", "sex"]]
df = pd.get_dummies(df, columns=["sequencing_batch", "rat_batch", "sex"], drop_first=True)

df.T.to_csv(outfile, sep="\t", header=False)
