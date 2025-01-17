import argparse
import pandas as pd
import numpy as np
import subprocess
import os

parser = argparse.ArgumentParser(
    description="Run GEMMA to get LOCO GRM PVE per gene"
)
parser.add_argument("pheno", help="BIMBAM phenotype file")
parser.add_argument("genes", help="gene IDs corresponding to phenotype file columns")
parser.add_argument("out", help="Output log file")
parser.add_argument("--chr", type=int, help="Chromosome (e.g. 1) to test")
parser.add_argument("--kinship", help="Kinship matrix file")
parser.add_argument("--covar", help="Covariates file")
parser.add_argument("--jobs", help="Number of gemma processes to run in parallel")
args = parser.parse_args()

WINDOW = 1e6

genes = pd.read_csv(args.genes, sep="\t")
genes = genes.to_dict("records")

run_genes = []
commands = []
for i, gene in enumerate(genes):
    if gene["chr"] != args.chr:
        continue
    command = [
        "gemma",
        "-p", f"<(cut -f{i + 1} {args.pheno})",
        "-vc", "1",
        "-k", args.kinship,
        "-outdir", f"{args.out}_tmp",
        "-o", gene['gene_id'],
    ]
    if args.covar is not None:
        command.extend(["-c", args.covar])
    command = " ".join(command)  # Get process substitution to work
    run_genes.append(gene["gene_id"])
    commands.append(command)
os.mkdir(f"{args.out}_tmp")
with open(f"{args.out}_tmp/commands.txt", "w") as f:
    for command in commands:
        f.write(command + "\n")
par_command = f"parallel -j {args.jobs} -a {args.out}_tmp/commands.txt"
subprocess.run(par_command, executable="/bin/bash", shell=True)

run_genes = [gene for gene in run_genes if os.path.exists(f"{args.out}_tmp/{gene}.log.txt")]
for i, gene in enumerate(run_genes):
    if i == 0:
        subprocess.run(f"cat {args.out}_tmp/{gene}.log.txt > {args.out}", shell=True)
    else:
        subprocess.run(f"cat {args.out}_tmp/{gene}.log.txt >> {args.out}", shell=True)
