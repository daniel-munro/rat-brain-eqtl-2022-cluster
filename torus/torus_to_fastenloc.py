import argparse
import pandas as pd
# import subprocess

p = argparse.ArgumentParser(description="Convert gene-variant pairs from tensorQTL to Matrix eQTL output format.")
p.add_argument("pipfile", help="Name of PIP file produced by torus")
p.add_argument("vcffile", help="Name of VCF file containing tested variants")
p.add_argument("tissue", help="Name of tissue (will be included in output).")
p.add_argument("outfile", help="Output file name (*.txt.gz)")
args = p.parse_args()

# pipfile = "~/br/data/dap/old/torus/pip.txt"
pip = pd.read_csv(args.pipfile, sep="\t", names = ["var", "gene_id", "PIP", "prior"])

# Sum PIP and num SNPs per "cluster" (gene):
pipclust = pip.groupby("gene_id", group_keys=False).apply(lambda x: x.PIP.sum())
pipclust = pipclust.reset_index(name="PIP_clust")
pipclust["PIP_clust"] = pipclust["PIP_clust"].clip(upper=1)
pip = pip.merge(pipclust, on="gene_id", how="left")
sizeclust = pip.groupby("gene_id", group_keys=False).apply(lambda x: x.shape[0])
sizeclust = sizeclust.reset_index(name="size_clust")
pip = pip.merge(sizeclust, on="gene_id", how="left")

pip = pip.loc[pip["PIP"] >= 1e-4, :] # Based on fastENLOC's summarize_dap2enloc.pl

snps = {}
for x in pip.to_dict(orient="records"):
    rec = f"[{x['gene_id']}:1@{args.tissue}={x['PIP']}[{x['PIP_clust']:g}:{x['size_clust']}]"
    if x["var"] in snps:
        snps[x["var"]].append(rec)
    else:
        snps[x["var"]] = [rec]

# vcffile = "~/br/data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
vcf = pd.read_csv(args.vcffile, sep="\t", comment="#", usecols=range(5), names=["chrom", "pos", "var", "ref", "alt"])
vcf = vcf.loc[vcf["var"].isin(snps)]
vcf["recs"] = ["|".join(snps[x]) for x in vcf["var"]]
vcf.to_csv(args.outfile, sep="\t", index=False, header=False)
