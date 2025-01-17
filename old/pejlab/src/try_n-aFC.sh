zcat ../data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz | grep -E '^#' > genotype.vcf
zcat ../data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz | grep -vE '^#' | sed 's/:/_/g' | sed -e 's/^/chr/g' >> genotype.vcf
bgzip genotype.vcf
tabix genotype.vcf.gz

zcat ../data/rsem_expected_count.gct.gz | tail -n +3 | cut -f1,3- | sed 's/\t/,/g' | sed 's/_Acbc//g' | gzip > rsem_n-aFC.txt.gz

python3 ../src/prepare_qtl_for_nafc.py ../data/fastqtl/main.genes.txt.gz > eqtls.txt

python3 n-aFC/src/nafc.py \
    --vcf genotype.vcf.gz \
    --expr rsem_n-aFC.txt.gz \
    --eqtl eqtls.txt \
    --output main.n-aFC.txt \
    --nthreads 16

python3 n-aFC/src/nafc.py \
    --vcf genotype.vcf.gz \
    --expr rsem_n-aFC.txt.gz \
    --eqtl eqtls.txt \
    --output conf.n-aFC.txt \
    --nthreads 16 \
    --conf

# subsets:
bcftools view -r chr12 genotype.vcf.gz -Oz -o genotype_subset.vcf.gz
tabix genotype_subset.vcf.gz
head -1 eqtls.txt > eqtls_subset.txt
grep "chr12_" eqtls.txt >> eqtls_subset.txt

python3 n-aFC/src/nafc.py \
    --vcf genotype_subset.vcf.gz \
    --expr rsem_n-aFC.txt.gz \
    --eqtl eqtls_subset.txt \
    --output main.n-aFC.txt

