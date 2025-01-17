import argparse
import pandas as pd
import numpy as np
import subprocess
from progress.bar import IncrementalBar
import os

parser = argparse.ArgumentParser(
    description="Run GEMMA on multiple genes with cis-windows"
)
parser.add_argument("geno", help="BIMBAM genotype file (*.gz)")
parser.add_argument("pheno", help="BIMBAM phenotype file")
parser.add_argument("snps", help="BIMBAM SNP annotation file")
parser.add_argument("genes", help="gene IDs corresponding to phenotype file columns")
parser.add_argument("out", help="Output file (TSV)")
parser.add_argument("--log", help="Output log file")
parser.add_argument("--chr", type=int, help="Chromosome (e.g. 1) to test")
parser.add_argument("--kinship", help="Kinship matrix file")
parser.add_argument("--covar", help="Covariates file")
args = parser.parse_args()

WINDOW = 1e6

# with open(args.genes, "r") as f:
#     genes = f.read().splitlines()
genes = pd.read_csv(args.genes, sep="\t")
genes = genes.loc[genes["chr"] == args.chr, :]
genes = genes.to_dict("records")

snps = pd.read_csv(args.snps, names=["snp", "pos", "chr"])  # , index_col="snp")
# snps = snps.sort_values(["chr", "pos"])
# snps = dict(tuple(snps.groupby("chr")))
# snps = snps.loc[snps["chr"] == 12, :]

bar = IncrementalBar("Genes tested", max=len(genes))
outfiles_created = False
for i, gene in enumerate(genes):
    # snp_first = snps["pos"].searchsorted(gene["pos"] - 1e6)
    # snp_last = snps["pos"].searchsorted(gene["pos"] + 1e6)
    in_window = np.where(
        np.logical_and(
            snps["chr"] == gene["chr"],
            np.logical_and(
                snps["pos"] >= gene["tss"] - WINDOW,
                snps["pos"] <= gene["tss"] + WINDOW,
            ),
        )
    )[0]
    if len(in_window) == 0:
        bar.next()
        continue
    first = in_window[0] + 1
    last = in_window[-1] + 1
    command = [
        "gemma",
        # process substitution doesn't work for genotype file:
        # "-g", f"<(head -{last} {args.geno} | tail -n +{first})",
        "-g", args.geno,
        "-snps", f"<(zcat {args.geno} | head -{last} | tail -n +{first} | cut -f1 -d',')",
        "-p", f"<(cut -f{i + 1} {args.pheno})",
        "-a", f"<(head -{last} {args.snps} | tail -n +{first})",
        "-outdir", f"{args.out}_tmp",
        "-o", gene['gene_id'],
        "-maf", "0.05",
    ]
    if args.covar is not None:
        command.extend(["-c", args.covar])
    if args.kinship is None:
        command.extend(["-lm", "1"])
    else:
        command.extend(["-lmm", "4", "-k", args.kinship])
    command = " ".join(command)  # Get process substitution to work
    # print(command)
    # stdout=subprocess.PIPE allows progress bar to print without GEMMA output.
    run = subprocess.run(command, executable="/bin/bash", shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Doesn't work:
    # if run.returncode != 0:
    #     print(f"Error running {gene['gene_id']}")

    # if os.path.isfile(f"output/{gene['gene_id']}.assoc.txt"):
    if os.path.isfile(f"{args.out}_tmp/{gene['gene_id']}.assoc.txt"):
        if not outfiles_created:
            subprocess.run(f"cat {args.out}_tmp/{gene['gene_id']}.assoc.txt | head -1 | sed 's/^/gene_id\t/' > {args.out}",
                        shell=True)
            subprocess.run(f"cat {args.out}_tmp/{gene['gene_id']}.log.txt > {args.log}",
                        shell=True)
            outfiles_created = True
        else:
            subprocess.run(f"cat {args.out}_tmp/{gene['gene_id']}.log.txt >> {args.log}",
                        shell=True)
        subprocess.run(f"cat {args.out}_tmp/{gene['gene_id']}.assoc.txt | tail -n +2 | sed 's/^/{gene['gene_id']}\t/' >> {args.out}",
                    shell=True)
    bar.next()
bar.finish()
