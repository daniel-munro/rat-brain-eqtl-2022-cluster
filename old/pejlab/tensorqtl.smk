rule vcf_to_plink:
    input:
        "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
    output:
        multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam")
    shell:
        """
        tools/plink2 --make-bed \
        --vcf {input} \
        --out "data/genotype/P50.rnaseq.88.unpruned"
        """

rule sample_participant_lookup:
    input: "rat_ids.txt"
    output: "data/expression/{region}.samples_participants.txt"
    run:
        with open(output[0], "w") as out:
            out.write("sample_id\tparticipant_id\n")
            for sample_id in sample_ids[wildcards.region]:
                rat_id = sample_id.split("_")[0]
                out.write("{0}\t{1}\n".format(sample_id, rat_id))

rule rsem_TPM:
    input:
        SFTP.remote("garibaldi.scripps.edu/gpfs/home/dmunro/data/rsem_out/{region}.rsem_TPM.gct.gz")
    output:
        "data/expression/{region}.rsem_TPM.gct.gz"
    shell:
        "cp {input} {output}"

rule rsem_expected_count:
    input:
        SFTP.remote("garibaldi.scripps.edu/gpfs/home/dmunro/data/rsem_out/{region}.rsem_expected_count.gct.gz")
    output:
        "data/expression/{region}.rsem_expected_count.gct.gz"
    shell:
        "cp {input} {output}"

rule prepare_expression:
    input:
        tpm_gct = "data/expression/{region}.rsem_TPM.gct.gz",
        counts_gct = "data/expression/{region}.rsem_expected_count.gct.gz",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
        samples = "data/expression/{region}.samples_participants.txt",
        chrlist = "data/genotype/P50.rnaseq.88.unpruned.chrlist.txt"
    output:
        "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
        "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi"
    params:
        prefix = "data/tensorqtl/{region}/{region}.main"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python3 src/eqtl_prepare_expression.py \
        {input.tpm_gct} {input.counts_gct} {input.gtf} \
        {input.samples} {input.chrlist} {params.prefix} \
        --tpm_threshold 0.1 \
        --count_threshold 6 \
        --sample_frac_threshold 0.2 \
        --normalization_method tmm
        """

rule expression_tsv_to_bed:
    input:
        tsv = "data/expression/ensemblGeneV2_Norm-counts-log2+1.txt",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        bed = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
        bedi = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi"
    params:
        interm = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed",
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        python3 src/tsv_to_bed.py {input.tsv} {input.gtf} {params.interm} --region {wildcards.region}
        bgzip {params.interm}
        tabix {output.bed}
        """

rule expression_gct_to_bed:
    input:
        gct = "data/expression/{region}.rsem_expected_count.gct.gz",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        bed = "data/expression/{region}.rsem_expected_count.bed.gz",
        bedi = "data/expression/{region}.rsem_expected_count.bed.gz.tbi"
    params:
        interm = "data/expression/{region}.rsem_expected_count.bed",
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        python3 src/gct_to_bed.py {input.gct} {input.gtf} {params.interm}
        bgzip {params.interm}
        tabix {output.bed}
        """

rule run_peer:
    input:
        # "data/tensorqtl/{region}/{region}.main.expression.bed.gz"
        "data/expression/{region}.rsem_expected_count.bed.gz",
    output:
        "data/tensorqtl/{region}/{region}.main.PEER_residuals.txt",
        "data/tensorqtl/{region}/{region}.main.PEER_alpha.txt",
        "data/tensorqtl/{region}/{region}.main.PEER_covariates.txt"
    params:
        prefix = "data/tensorqtl/{region}/{region}.main",
        num_peer = 15
    conda:
        "envs/peer.yaml"
    shell:
        "Rscript src/run_PEER.R {input} {params.prefix} {params.num_peer}"

rule run_peer_main2:
    input:
        "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz"
    output:
        "data/tensorqtl/{region}/{region}.main2.PEER_residuals.txt",
        "data/tensorqtl/{region}/{region}.main2.PEER_alpha.txt",
        "data/tensorqtl/{region}/{region}.main2.PEER_covariates.txt"
    params:
        prefix = "data/tensorqtl/{region}/{region}.main2",
        num_peer = 15
    conda:
        "envs/peer.yaml"
    shell:
        "Rscript src/run_PEER.R {input} {params.prefix} {params.num_peer}"

rule other_covariates:
    input:
        "data/DemographicInfo_TN_RNASeq_88_LauraSaba.csv"
    output:
        "data/tensorqtl/other_covariates.tsv"
    shell:
        "python3 src/tensorqtl_covar.py {input} {output}"

rule similarity_to_founders:
    input:
        pop_vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        pop_vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        founder_vcf = "data/genotype/founders.vcf.gz",
        founder_vcfi = "data/genotype/founders.vcf.gz.tbi"
    output:
        "data/tensorqtl/sim_to_founders.txt"
    conda:
        "envs/biopython.yaml"
    shell:
        "python3 src/sim_to_each_founder.py {input.pop_vcf} {input.founder_vcf} {output}"

rule combine_covariates:
    input:
        peer = "data/tensorqtl/{region}/{region}.{method}.PEER_covariates.txt",
        founder_sim = "data/tensorqtl/sim_to_founders.txt",
        other = "data/tensorqtl/other_covariates.tsv"
    output:
        "data/tensorqtl/{region}/{region}.{method}.combined_covariates.txt"
    params:
        prefix = "data/tensorqtl/{region}/{region}.{method}"
    shell:
        """
        python3 src/combine_covariates.py {input.peer} {params.prefix} \
        --add_covariates {input.founder_sim} {input.other}
        """

def expr_bed(wildcards):
    return {
        # "basic": "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
        "basic": "data/expression/{region}.rsem_expected_count.bed.gz",
        "basic2": "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz"
    }[wildcards.method]

rule empty_covariates:
    """tensorQTL currently requires a covariates file."""
    input:
        # "data/expression/ensemblGeneV2_Norm-counts-log2+1.bed.gz"
        expr_bed
    output:
        "data/tensorqtl/{region}/{region}.{method}.covar_empty.txt"
    shell:
        "zcat {input} | head -1 | sed 's/gene_id/ID/' | cut -f4- > {output} || true"

rule run_tensorqtl_perm_main:
    input:
        geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
        # bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
        # bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
        bed = "data/expression/{region}.rsem_expected_count.bed.gz",
        bedi = "data/expression/{region}.rsem_expected_count.bed.gz.tbi",
        covar = "data/tensorqtl/{region}/{region}.main.combined_covariates.txt"
    output:
        "data/tensorqtl/{region}/{region}.main.cis_qtl.txt.gz"
    conda:
        "envs/tensorqtl.yaml"
    shell:
        """
        # pip install -e tools/tensorqtl
        python -m tensorqtl --mode cis \
        data/genotype/P50.rnaseq.88.unpruned \
        {input.bed} \
        {wildcards.region}.main \
        --covariates {input.covar} \
        --output_dir data/tensorqtl/{wildcards.region}
        """

rule run_tensorqtl_perm_main2:
    input:
        geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
        bed = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
        bedi = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi",
        covar = "data/tensorqtl/{region}/{region}.main2.combined_covariates.txt"
    output:
        "data/tensorqtl/{region}/{region}.main2.cis_qtl.txt.gz"
    conda:
        "envs/tensorqtl.yaml"
    shell:
        """
        # pip install -e tools/tensorqtl
        python -m tensorqtl --mode cis \
        data/genotype/50.rnaseq.88.unpruned \
        {input.bed} \
        {wildcards.region}.main2 \
        --covariates {input.covar} \
        --output_dir data/tensorqtl/{wildcards.region}
        """

rule run_tensorqtl_perm_basic:
    input:
        geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
        # bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
        # bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
        bed = "data/expression/{region}.rsem_expected_count.bed.gz",
        bedi = "data/expression/{region}.rsem_expected_count.bed.gz.tbi",
        covar = "data/tensorqtl/{region}/{region}.basic.covar_empty.txt"
    output:
        "data/tensorqtl/{region}/{region}.basic.cis_qtl.txt.gz"
    conda:
        "envs/tensorqtl.yaml"
    shell:
        """
        # pip install -e tools/tensorqtl
        python -m tensorqtl --mode cis \
        data/genotype/P50.rnaseq.88.unpruned \
        {input.bed} \
        {wildcards.region}.basic \
        --covariates {input.covar} \
        --output_dir data/tensorqtl/{wildcards.region}
        """

rule run_tensorqtl_perm_basic2:
    input:
        geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
        bed = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
        bedi = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi",
        covar = "data/tensorqtl/{region}/{region}.basic2.covar_empty.txt"
    output:
        "data/tensorqtl/{region}/{region}.basic2.cis_qtl.txt.gz"
    conda:
        "envs/tensorqtl.yaml"
    shell:
        """
        # pip install -e tools/tensorqtl
        python -m tensorqtl --mode cis \
        data/genotype/P50.rnaseq.88.unpruned \
        {input.bed} \
        {wildcards.region}.basic2 \
        --covariates {input.covar} \
        --output_dir data/tensorqtl/{wildcards.region}
        """
