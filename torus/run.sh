# TISSUE=Acbc
# TENSORPREFIX=AQCT
TISSUE=$1
TENSORPREFIX=$2

# python tensorqtl_to_matrixeqtl.py \
#     ../tensorqtl/$TENSORPREFIX.cis_qtl.txt.gz \
#     ../tensorqtl/$TENSORPREFIX/ \
#     mateqtl.$TISSUE.txt.gz

# torus/src/torus -d mateqtl.$TISSUE.txt.gz -dump_pip pip.$TISSUE.txt

python torus_to_fastenloc.py \
    pip.$TISSUE.txt \
    ../genotype/P50.rnaseq.88.unpruned.vcf.gz \
    $TISSUE \
    fastenloc.eqtl.anno.$TISSUE.vcf

bgzip fastenloc.eqtl.anno.$TISSUE.vcf
