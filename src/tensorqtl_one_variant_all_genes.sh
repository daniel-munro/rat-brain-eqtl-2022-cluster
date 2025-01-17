#!/usr/bin/env bash
# set -euo pipefail

# TensorQTL has trans mode, but it filters out variants in cis region, always +/- 5Mbp.
# So do that plus cis_nominal mode with 5Mbp window.

VAR="chr16:76543372"
VAR2="chr16_76543372"
CHR=16
BED=eyes/data/tensorqtl/eyes.expression.bed.gz
NEWBED=eyes/data/tensorqtl/one_variant/${CHR}.bed
PREFIX=$VAR2
COVAR=eyes/data/tensorqtl/main4.PEER_covariates.txt
GENO=eyes/data/genotype/eyes
NEWGENO=eyes/data/tensorqtl/one_variant/$VAR2
OUTDIR=eyes/data/tensorqtl/one_variant

# # Necessary because there's an error if there are no variants on a gene's chromosome in cis_nominal mode.
# zcat $BED | head -1 > $NEWBED
# zcat $BED | grep -P "^${CHR}\t" >> $NEWBED
# bgzip $NEWBED

# plink2 --bfile $GENO --extract <(echo $VAR) --make-bed --out $NEWGENO

python -m tensorqtl --mode cis_nominal \
    $NEWGENO \
    $NEWBED.gz \
    $PREFIX \
    --covariates $COVAR \
    --window 5000000 \
    --output_dir $OUTDIR

# --return_dense doesn't work for some reason.
python -m tensorqtl --mode trans \
    $NEWGENO \
    $BED \
    $PREFIX \
    --covariates $COVAR \
    --pval_threshold 1.00001 \
    --output_dir $OUTDIR

    #--pval_threshold 1.0 \
# Then use pandas to convert parquet to tsv, but it isn't working.
