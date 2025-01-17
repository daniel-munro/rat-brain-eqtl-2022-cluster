import argparse
import pandas as pd
import numpy as np
import subprocess
import os

def kinship_command(gene, snps, args):
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
        return None
    first = in_window[0] + 1
    last = in_window[-1] + 1
    snps_inp = f"<(zcat {args.geno} | head -{last} | tail -n +{first} | cut -f1 -d',')"
    a = f"<(head -{last} {args.snps} | tail -n +{first})"
    outdir = f"{args.tmpdir}/kinship"
    out = gene['gene_id']
    command = f"gemma -gk -g {args.geno} -snps {snps_inp} -a {a} -p {args.pheno} -outdir {outdir} -o {out}"
    return command

def assoc_command(i, gene, args):
    pheno = f"<(cut -f{i + 1} {args.pheno})"
    kinship = f"{args.tmpdir}/kinship/{gene['gene_id']}.cXX.txt"
    outdir = f"{args.tmpdir}/assoc"
    out = gene['gene_id']
    command = f"gemma -vc 1 -p {pheno} -k {kinship} -outdir {outdir} -o {out}"
    if args.covar is not None:
        command += f" -c {args.covar}"
    return command

parser = argparse.ArgumentParser(
    description="Run GEMMA on to get cis-window PVE per gene"
)
parser.add_argument("geno", help="BIMBAM genotype file (*.gz)")
parser.add_argument("pheno", help="BIMBAM phenotype file")
parser.add_argument("snps", help="BIMBAM SNP annotation file")
parser.add_argument("genes", help="gene IDs corresponding to phenotype file columns")
parser.add_argument("out", help="Output file (TSV)")
parser.add_argument("tmpdir", help="Directory for temp files")
parser.add_argument("--covar", help="Covariates file")
parser.add_argument("--chr", type=int, help="Chromosome (e.g. 1) to test")
parser.add_argument("--jobs", type=int, default=16, help="Number of gemma processes to run in parallel")
args = parser.parse_args()

WINDOW = 1e6

genes = pd.read_csv(args.genes, sep="\t")
genes = genes.to_dict("records")

snps = pd.read_csv(args.snps, names=["snp", "pos", "chr"])

run_genes = []
commands = []
for i, gene in enumerate(genes):
    if gene["chr"] != args.chr:
        continue
    command1 = kinship_command(gene, snps, args)
    if command1 is None:
        continue
    command2 = assoc_command(i, gene, args)
    command = f"{command1} && {command2}"
    run_genes.append(gene["gene_id"])
    commands.append(command)

os.mkdir(f"{args.tmpdir}/kinship")
os.mkdir(f"{args.tmpdir}/assoc")
with open(f"{args.tmpdir}/commands.txt", "w") as f:
   for command in commands:
       f.write(command + "\n")
# par_command = f"parallel -j {args.jobs} -a {args.tmpdir}/commands.txt"
# For some reason the above didn't work but this does!:
par_command = f"parallel --version && parallel -j {args.jobs} -a {args.tmpdir}/commands.txt"
# par_command = f"bash {args.tmpdir}/commands.txt"
subprocess.run(par_command, executable="/usr/bin/bash", shell=True)
# par_command = f"parallel -j {args.jobs} -a {args.tmpdir}/commands.txt"
# subprocess.run(par_command, executable="/usr/bin/bash", shell=True)

run_genes = [gene for gene in run_genes if os.path.exists(f"{args.tmpdir}/assoc/{gene}.log.txt")]
for i, gene in enumerate(run_genes):
    redir = ">" if i == 0 else ">>"
    subprocess.run(f"cat {args.tmpdir}/assoc/{gene}.log.txt {redir} {args.out}", shell=True)
