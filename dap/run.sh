DSET=$1

# python3 prepare_sbams.py \
#         ../genotype/P50.rnaseq.88.unpruned.vcf.gz \
#         ../expression/ensembl-gene_inv-quant_$DSET.bed.gz \
#         input/$DSET/ \
#         --covar covar/$DSET.covariates.txt

# zcat ../expression/ensembl-gene_inv-quant_$DSET.bed.gz | grep -v "^#" | cut -f4 > genes.$DSET.txt

parallel -j8 sh run_dap_one_gene.sh ::: $DSET :::: genes.$DSET.txt

# perl ~/tools/fastenloc/src/summarize_dap2enloc.pl \
#         -dir output/$DSET/ \
#         -vcf ../genotype/P50.rnaseq.88.unpruned.vcf.gz \
#         -tissue $DSET \
#         | bgzip > fastenloc.eqtl.annotation.$DSET.vcf.gz
