import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description="Shuffle expression files to mix brain regions")
parser.add_argument("-i", nargs="+", help="Expression bed files for all regions")
parser.add_argument("-o", nargs="+", help="Output bed files, will have same sample counts as inputs")
parser.add_argument("--seed", type=int, help="random number generator seed so output is replicable")
args = parser.parse_args()

beds = [pd.read_csv(f, sep="\t", index_col=3, dtype=str) for f in args.i]
# Take intersection of genes from all regions and use it for all outputs.
# (The original expression txt files also have different gene sets FYI).
genes = set.intersection(*[set(bed.index) for bed in beds])
genes = [gene for gene in beds[0].index if gene in genes] # Keep original gene order
beds = [bed.loc[genes, :] for bed in beds]

indivs = sorted(set.union(*[set(bed.columns[3:]) for bed in beds]))
insamples = {ind: [] for ind in indivs}
for region, bed in enumerate(beds):
    for ind in bed.columns[3:]:
        insamples[ind].append(region)

outbeds = [bed.iloc[:, :3].copy() for bed in beds]
rng = np.random.default_rng(seed=args.seed)
for ind in indivs:
    outsamples = rng.permutation(insamples[ind])
    for i, region in enumerate(insamples[ind]):
        outbeds[outsamples[i]][ind] = beds[region][ind]

for i, outbed in enumerate(outbeds):
    outbed.insert(3, "gene_id", outbed.index)
    outbed.to_csv(args.o[i], sep="\t", index=False)
