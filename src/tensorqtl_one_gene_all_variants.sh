#!/usr/bin/env bash
# set -euo pipefail

## I copied this to the eyes directory.
# # IOP (eye) QTL:
# GENE=ENSRNOG00000016696
# BED=eyes/data/tensorqtl/eyes.expression.bed.gz
# NEWBED=eyes/data/tensorqtl/one_gene/$GENE.bed
# PREFIX=$GENE
# COVAR=eyes/data/tensorqtl/main4.PEER_covariates.txt
# GENO=eyes/data/genotype/eyes
# OUTDIR=eyes/data/tensorqtl/one_gene
# 
# # zcat $BED | head -1 > $NEWBED
# # zcat $BED | grep $GENE >> $NEWBED
# # bgzip $NEWBED
# 
# # python -m tensorqtl --mode cis_nominal \
# #     $GENO \
# #     $NEWBED.gz \
# #     $PREFIX \
# #     --covariates $COVAR \
# #     --output_dir $OUTDIR
# 
# # Larger window for figure:
# python -m tensorqtl --mode cis_nominal \
#     $GENO \
#     $NEWBED.gz \
#     $PREFIX.4mbp \
#     --covariates $COVAR \
#     --window 4000000 \
#     --output_dir $OUTDIR


# # one-kidney QTL for Joel:
# GENE=ENSRNOG00000002227
# BED=data/expression/ensembl-gene_inv-quant_ComBat_Acbc.bed.gz
# NEWBED=data/tensorqtl/one_gene/$GENE.bed
# PREFIX=$GENE
# COVAR=data/tensorqtl/AECT.combat_covariates.txt
# GENO=data/genotype/P50.rnaseq.88.unpruned
# OUTDIR=data/tensorqtl/one_gene

# zcat $BED | head -1 > $NEWBED
# zcat $BED | grep $GENE >> $NEWBED
# bgzip $NEWBED
# # tabix $NEWBED

# # Larger window to include SNPs in question:
# python -m tensorqtl --mode cis_nominal \
#     $GENO \
#     $NEWBED.gz \
#     $PREFIX \
#     --covariates $COVAR \
#     --window 1500000 \
#     --output_dir $OUTDIR

# Then use pandas to convert parquet to tsv, but it isn't working.
