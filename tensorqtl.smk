localrules:
    plink_sample_list,
    vcf_to_plink,
    expression_tsv_to_bed,
    expression_pcs,
    combine_covariates,
    empty_covariates,


rule plink_sample_list:
    """Assumes any expression file has same individuals as this one for that tissue."""
    input:
        "data/expression/{base}.bed.gz"
    output:
        "data/tensorqtl/genotypes/{base}.samples.txt"
    shell:
        "zcat {input} | head -1 | cut -f5- | sed 's/\\t/\\n/g' > {output} || true"


def geno_prefix(wildcards):
    """Genotypes are expression dataset-specific to remove monomorphic SNPs."""
    tissue = dict(I="IL", L="LHb", N="Acbc", O="VoLo", P="PL")[wildcards.code[0]]
    expr = dict(L="log2", Q="inv-quant")[wildcards.code[1]]
    return f"data/tensorqtl/genotypes/ensembl-gene_{expr}_{tissue}"


def expr_bed(wildcards):
    tissue = dict(I="IL", L="LHB", N="Acbc", O="VoLo", P="PL")[wildcards.code[0]]
    expr = dict(L="log2", Q="inv-quant")[wildcards.code[1]]
    return f"data/expression/ensembl-gene_{expr}_{tissue}.bed.gz"


def expr_pca_file(wildcards):
    tissue = dict(I="IL", L="LHB", N="Acbc", O="VoLo", P="PL")[wildcards.code[0]]
    expr = dict(L="log2", Q="inv-quant")[wildcards.code[1]]
    return f"data/expression/pca/ensembl-gene_{expr}_PCs_{tissue}.20.txt"


def covar_file(wildcards):
    covar = dict(N="covar_empty", C="combined_covariates")[wildcards.code[2]]
    return f"data/tensorqtl/{wildcards.code}.{covar}.txt"


rule vcf_to_plink:
    """Get SNPs that are not monomorphic in a given set of samples."""
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        samples = "data/tensorqtl/genotypes/{base}.samples.txt"
    output:
        multiext("data/tensorqtl/genotypes/{base}", ".bed", ".bim", ".fam")
    params:
        prefix = "data/tensorqtl/genotypes/{base}"
    shell:
        # Use these for all monomorphic SNPs:
        # --maf 0.0001 \
        # --max-maf 0.9999 \
        """
        plink2 --make-bed \
        --vcf {input.vcf} \
        --keep {input.samples} \
        --maf 0.01 \
        --max-maf 0.99 \
        --out {params.prefix}
        """


rule expression_tsv_to_bed:
    input:
        tsv = "data/expression/ensembl-gene_{expr}_{tissue}.txt",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        bed = "data/expression/ensembl-gene_{expr}_{tissue}.bed.gz",
        bedi = "data/expression/ensembl-gene_{expr}_{tissue}.bed.gz.tbi"
    params:
        interm = "data/expression/ensembl-gene_{expr}_{tissue}.bed"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        """
        python3 src/tsv_to_bed.py {input.tsv} {input.gtf} {params.interm}
        bgzip {params.interm}
        tabix {output.bed}
        """


rule expression_pcs:
    input:
        "data/expression/{base}.txt"
    output:
        "data/expression/{base}.20.txt"
    run:
        d = pd.read_csv(input[0], sep="\t", dtype="str")
        d.index = [x.split("_")[0] for x in d.index.to_list()]
        d = d.transpose().iloc[:20, :]
        d.to_csv(output[0], sep="\t", index_label="ID")


rule similarity_to_founders:
    input:
        pop_vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        pop_vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        founder_vcf = "data/genotype/founders.vcf.gz",
        founder_vcfi = "data/genotype/founders.vcf.gz.tbi"
    output:
        "data/tensorqtl/sim_to_founders.txt"
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python3 src/sim_to_each_founder.py {input.pop_vcf} {input.founder_vcf} {output}
        sed -i 's/-N//g' {output}
        """


rule combine_covariates:
    input:
        pca = expr_pca_file,
        founder_sim = "data/tensorqtl/sim_to_founders.txt",
    output:
        "data/tensorqtl/{code}.combined_covariates.txt"
    params:
        prefix = "data/tensorqtl/{code}"
    shell:
        """
        python3 src/combine_covariates.py {input.pca} {params.prefix} \
        --add_covariates {input.founder_sim}
        """


rule empty_covariates:
    """tensorQTL currently requires a covariates file."""
    input:
        expr_bed
    output:
        "data/tensorqtl/{code}.covar_empty.txt"
    shell:
        "zcat {input} | head -1 | sed 's/gene_id/ID/' | cut -f4- > {output} || true"


rule tensorqtl_perm:
    input:
        geno = lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam"),
        bed = expr_bed,
        bedi = lambda w: expr_bed(w) + ".tbi",
        covar = covar_file
    output:
        "data/tensorqtl/{code}.cis_qtl.txt.gz"
    params:
        geno_prefix = geno_prefix,
        outdir = "data/tensorqtl"
    # conda:
    #     "../envs/tensorqtl.yaml"
    resources:
        walltime = 12,
        partition = "--partition=gpu",
    shell:
        """
        module load cuda
        python -m tensorqtl \
            --mode cis \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.code} \
            --covariates {input.covar} \
            --output_dir {params.outdir}
        """


rule tensorqtl_independent:
    input:
        geno = lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam"),
        bed = expr_bed,
        bedi = lambda w: expr_bed(w) + ".tbi",
        covar = covar_file,
        cis = "data/tensorqtl/{code}.cis_qtl.txt.gz"
    output:
        "data/tensorqtl/{code}.cis_independent_qtl.txt.gz"
    params:
        geno_prefix = geno_prefix,
        outdir = "data/tensorqtl"
    # conda:
    #     "../envs/tensorqtl.yaml"
    resources:
        walltime = 12,
        partition = "--partition=gpu",
    shell:
        """
        module load cuda
        python -m tensorqtl \
            --mode cis_independent \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.code} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            --output_dir {params.outdir}
        """


rule tensorqtl_nominal:
    input:
        geno = lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam"),
        bed = expr_bed,
        bedi = lambda w: expr_bed(w) + ".tbi",
        covar = covar_file
    output:
        expand("data/tensorqtl/{{code}}/{{code}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    params:
        geno_prefix = geno_prefix,
        outdir = "data/tensorqtl/{code}"
    # conda:
    #     "../envs/tensorqtl.yaml"
    resources:
        walltime = 12,
        partition = "--partition=gpu",
    shell:
        """
        # module load cuda
        # pip install -e ~/tools/tensorqtl
        mkdir -p {params.outdir}
        python -m tensorqtl --mode cis_nominal \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.code} \
            --covariates {input.covar} \
            --output_dir {params.outdir}
        """


rule tensorqtl_all_signif:
    input:
        perm = "data/tensorqtl/{code}.cis_qtl.txt.gz",
        nom = expand("data/tensorqtl/{{code}}/{{code}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    output:
        "data/tensorqtl/{code}.cis_qtl_signif.txt.gz"
    params:
        nom_prefix = "data/tensorqtl/{code}/{code}"
    # conda:
    #     "../envs/tensorqtl.yaml"
    shell:
        "python3 src/tensorqtl_all_sig_eqtls.py {input.perm} {params.nom_prefix} {output}"


rule tensorqtl_trans:
    input:
        geno = lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam"),
        bed = expr_bed,
        bedi = lambda w: expr_bed(w) + ".tbi",
        covar = covar_file
    output:
        "data/tensorqtl/{code}.trans_qtl_pairs.txt.gz"
    params:
        geno_prefix = "data/genotype/P50.rnaseq.88.unpruned",
    # conda:
        # "../envs/tensorqtl.yaml"
    resources:
        walltime = 12,
        # partition = "--partition=gpu",
    run:
        # Must be run from tensorqtl environment.
        # Run in python script so cis variants aren't filtered out.
        # Based on examples at https://github.com/broadinstitute/tensorqtl
        import tensorqtl
        from tensorqtl import genotypeio, trans
        pheno, pheno_pos = tensorqtl.read_phenotype_bed(input.bed)
        covar = pd.read_csv(input.covar, sep="\t", index_col=0).T
        geno = genotypeio.PlinkReader(params.geno_prefix).load_genotypes()
        # Doesn't work due to bug (get_all_genotypes not defined):
        # geno = tensorqtl.read_genotypes(params.geno_prefix)
        d = trans.map_trans(geno, pheno, covar)
        d.to_csv(output[0], sep="\t", index=False, float_format="%.6g")


# rule tensorqtl_perm_chr:
#     """For debugging."""
#     input:
#         # geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
#         geno = lambda w: multiext(geno_prefix(w), ".bed", ".bim", ".fam"),
#         bed = expr_bed,
#         bedi = lambda w: expr_bed(w) + ".tbi",
#         covar = covar_file
#     output:
#         "data/tensorqtl/{code}.{chrn}.cis_qtl.txt.gz"
#     params:
#         geno_prefix = "data/genotype/P50.rnaseq.88.unpruned",
#         outdir = "data/tensorqtl"
#     conda:
#         "envs/tensorqtl.yaml"
#     resources:
#         walltime = 12,
#     #     # gpu = ",nodes=1:gpus=1"
#         # gpu = ",nodes=1:ppn=8:gtx1080ti -q gpu"
#     #     # gpu = ",nodes=1:ppn=20:gtx1080 -q gpu"
#     shell:
#         """
#         zcat {input.bed} | head -1 > {params.outdir}/{code}.{chrn}.bed
#         zcat {input.bed} | ......
#         python -m tensorqtl --mode cis \
#             {params.geno_prefix} \
#             {input.bed} \
#             {wildcards.code}.{wildcards.chrn} \
#             --covariates {input.covar} \
#             --output_dir {params.outdir}
#         """
