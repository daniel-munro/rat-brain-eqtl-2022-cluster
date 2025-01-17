import argparse
import pandas as pd
import numpy as np
import subprocess
# from progress.bar import IncrementalBar
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
parser.add_argument("--cis_only", action="store_true", help="If provided, only test cis-window (+/- 1Mbp) per gene")
parser.add_argument("--nsnps", type=int, required=False, help="Test only random subset of N SNPs per gene (different set per gene). E.g. for QQ plots. Incompatible with cis_only.")
parser.add_argument("--pval_threshold", type=float, required=False, help="Maximum p-value to include in combined output")
parser.add_argument("--jobs", type=int, default=16, help="Number of gemma processes to run in parallel")
parser.add_argument("--continue", action="store_true", dest="cont", help="If provided, genes with existing output files will not be rerun.")
args = parser.parse_args()

WINDOW = 1e6

genes = pd.read_csv(args.genes, sep="\t")
# genes = genes.loc[genes["chr"] == args.chr, :] NO, keep original list to correspond to pheno file columns.
genes = genes.to_dict("records")

snps = pd.read_csv(args.snps, names=["snp", "pos", "chr"])  # , index_col="snp")
# snps = snps.sort_values(["chr", "pos"])
# snps = dict(tuple(snps.groupby("chr")))
# snps = snps.loc[snps["chr"] == 12, :]

# bar = IncrementalBar("Genes tested", max=len(genes))
# outfiles_created = False
run_genes = []
commands = []

for i, gene in enumerate(genes):
    if gene["chr"] != args.chr:
        continue
    if args.cis_only:
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
            # bar.next()
            continue
        first = in_window[0] + 1
        last = in_window[-1] + 1
        snps_input = f"<(zcat {args.geno} | head -{last} | tail -n +{first} | cut -f1 -d',')"
        a_input = f"<(head -{last} {args.snps} | tail -n +{first})"
    elif args.nsnps is not None:
        lines = np.sort(np.random.choice(snps.shape[0], size=args.nsnps, replace=False) + 1)
        sed_command = ";".join([f"{x}p" for x in lines])
        snps_input = f"<(zcat {args.geno} | sed -n '{sed_command}' | cut -f1 -d',')"
        a_input = f"<(sed -n '{sed_command}' {args.snps})"
    else:
        snps_input = f"<(zcat {args.geno} | cut -f1 -d',')"
        a_input = args.snps

    command = [
        "gemma",
        # process substitution doesn't work for genotype file:
        # "-g", f"<(head -{last} {args.geno} | tail -n +{first})",
        "-g", args.geno,
        "-snps", snps_input,
        "-p", f"<(cut -f{i + 1} {args.pheno})",
        "-a", a_input,
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
    run_genes.append(gene["gene_id"])
    if not (args.cont and os.path.exists(f"{args.out}_tmp/{gene['gene_id']}.assoc.txt")):
        commands.append(command)

# If not in continue mode, throw an error if directory exists so files aren't overwritten.
if (not args.cont) or (not os.path.exists(f"{args.out}_tmp")):
    os.mkdir(f"{args.out}_tmp")
with open(f"{args.out}_tmp/commands.txt", "w") as f:
   for command in commands:
       f.write(command + "\n")
par_command = f"parallel -j {args.jobs} -a {args.out}_tmp/commands.txt"
subprocess.run(par_command, executable="/bin/bash", shell=True)

run_genes = [gene for gene in run_genes if os.path.exists(f"{args.out}_tmp/{gene}.assoc.txt")]
for i, gene in enumerate(run_genes):
    if i == 0:
        subprocess.run(f"cat {args.out}_tmp/{gene}.assoc.txt | head -1 | sed 's/^/gene_id\t/' > {args.out}",
                    shell=True)
        subprocess.run(f"cat {args.out}_tmp/{gene}.log.txt > {args.log}",
                    shell=True)
    else:
        subprocess.run(f"cat {args.out}_tmp/{gene}.log.txt >> {args.log}",
                    shell=True)
    pval_col = 11 if args.kinship is None else 13
    pfilter = "" if args.pval_threshold is None else f"awk '${pval_col} < {args.pval_threshold}' |"
    subprocess.run(f"cat {args.out}_tmp/{gene}.assoc.txt | tail -n +2 | {pfilter} sed 's/^/{gene}\t/' >> {args.out}",
                shell=True)
    # bar.next()
# bar.finish()
