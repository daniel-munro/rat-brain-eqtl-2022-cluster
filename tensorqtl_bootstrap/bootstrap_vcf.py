import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Make bootstrap sampling of VCF file, duplicating samples where necessary."
)
parser.add_argument("vcf", help="VCF file")
parser.add_argument("samples", help="File containing sampled list of samples")
parser.add_argument("out", help="Output file (VCF but with no header except column names)")
args = parser.parse_args()

# skiprows is specific to my file. I can't use comment = '##' because only 1 char allowed, and column header row starts with '#'.
vcf = pd.read_csv(args.vcf, sep="\t", skiprows=34)
samples = list(pd.read_csv(args.samples, names=["sample"])["sample"])

vcf2 = vcf.iloc[:, :9].copy()
count = {sample: 0 for sample in samples}
for sample in samples:
    count[sample] += 1
    vcf2[f"{sample}_{count[sample]}"] = vcf[sample]

vcf2.to_csv(args.out, sep="\t", index=False)
