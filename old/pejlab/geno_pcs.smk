rule prepare_smartpca_vcf_to_plink:
    # Convert to PLINK2 files.
    input:
        "data/genotype/round8_unpruned.vcf.gz"
    output:
        expand("data/genotype_pcs/round8_unpruned.{ext}", ext=["bed", "bim", "fam", "log"])
    params:
        prefix1 = "data/genotype_pcs/round8_unpruned"
    shell:
        # See https://github.com/broadinstitute/gtex-pipeline/blob/master/genotype/compute_genotype_pcs.py
        """
        tools/plink2 --vcf {input} \
        --out {params.prefix1} \
        --make-bed --max-alleles 2 \
        --output-chr chrM
        """

rule prepare_smartpca_filter_maf05_geno01:
    # Filter by minor allele frequency and missingness.
    input:
        expand("data/genotype_pcs/round8_unpruned.{ext}", ext=["bed", "bim", "fam"])
    output:
        expand("data/genotype_pcs/round8_unpruned.maf05_geno01.{ext}", ext=["bed", "bim", "fam", "log"])
    params:
        prefix1 = "data/genotype_pcs/round8_unpruned",
        prefix2 = "data/genotype_pcs/round8_unpruned.maf05_geno01"
    shell:
        """
        tools/plink2 --make-bed --output-chr chrM \
        --bfile {params.prefix1} \
        --maf 0.05 --geno 0.01 \
        --out {params.prefix2}
        """

rule prepare_smartpca_prune_LD:
    # Get independent SNPs.
    input:
        expand("data/genotype_pcs/round8_unpruned.maf05_geno01.{ext}", ext=["bed", "bim", "fam"])
    output:
        expand("data/genotype_pcs/round8_unpruned.maf05_geno01.prune.{ext}", ext=["in", "out"])
    params:
        prefix2 = "data/genotype_pcs/round8_unpruned.maf05_geno01"
    shell:
        """
        tools/plink2 --bfile {params.prefix2} \
        --indep-pairwise 200 100 0.1 \
        --out {params.prefix2}
        """

rule prepare_smartpca_subset_to_pruned:
    # Subset to independent SNPs. Output to pgen first since --sort-vars is not yet supported for bed in PLINK 2.
    input:
        # "data/genotype_pcs/round8_unpruned.maf05_geno01.bed"
        expand("data/genotype_pcs/round8_unpruned.maf05_geno01.{ext}", ext=["bed", "bim", "fam", "prune.in"])
    output:
        expand("data/genotype_pcs/round8_unpruned.maf05_geno01.pruned.{ext}", ext=["pgen", "psam", "pvar", "log"])
    params:
        prefix2 = "data/genotype_pcs/round8_unpruned.maf05_geno01",
        prefix3 = "data/genotype_pcs/round8_unpruned.maf05_geno01.pruned"
    shell:
        """
        tools/plink2 --bfile {params.prefix2} \
        --output-chr chrM \
        --extract {params.prefix2}.prune.in \
        --out {params.prefix3} \
        --sort-vars \
        --make-pgen        
        """

rule prepare_smartpca_pruned_to_bed:
    # Convert pgen/psam/pvar to bed/bim/fam for smartpca.
    input:
        expand("data/genotype_pcs/round8_unpruned.maf05_geno01.pruned.{ext}", ext=["pgen", "psam", "pvar"])
    output:
        bed = "data/genotype_pcs/round8_unpruned.maf05_geno01.pruned.bed",
        bim = "data/genotype_pcs/round8_unpruned.maf05_geno01.pruned.bim",
        fam = "data/genotype_pcs/round8_unpruned.maf05_geno01.pruned.fam"
    params:
        prefix3 = "data/genotype_pcs/round8_unpruned.maf05_geno01.pruned"
    shell:
        """
        tools/plink2 --pfile {params.prefix3} \
        --output-chr chrM \
        --make-bed \
        --out {params.prefix3}
        sed --in-place='' 's/-9/1/g' {output.fam}
        """

rule run_smartpca:
    input:
        bed = "data/genotype_pcs/round8_unpruned.maf05_geno01.pruned.bed",
        bim = "data/genotype_pcs/round8_unpruned.maf05_geno01.pruned.bim",
        fam = "data/genotype_pcs/round8_unpruned.maf05_geno01.pruned.fam"
    output:
        pca = "data/genotype_pcs/smartpca.pca",
        evl = "data/genotype_pcs/smartpca.eval",
        evec = "data/genotype_pcs/smartpca.pca.evec",
        log = "data/genotype_pcs/smartpca.log"
    params:
        plot = "data/genotype_pcs/smartpca.plot" # Not produced, I think since eigplot program not found.
    # conda:
    #     "envs/smartpca.yaml"
    shell:
        """
        PATH=tools/EIG-7.2.1/bin:$PATH smartpca.perl \
        -i {input.bed} \
        -a {input.bim} \
        -b {input.fam} \
        -k 20 -m 0 \
        -o {output.pca} \
        -e {output.evl} \
        -p {params.plot} \
        -l {output.log}
        """

rule prepare_genotype_pcs:
    input:
        "data/genotype_pcs/smartpca.pca.evec"
    output:
        "data/genotype_pcs/genotype_pcs.tsv"
    params:
        n_pcs = 3
    shell:
        "python3 src/prepare_genotype_pcs.py {input} {output} {params.n_pcs}"
