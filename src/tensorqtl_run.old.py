# Run in python script to use random_tiebreak option
# Based on examples at https://github.com/broadinstitute/tensorqtl
import argparse
import pandas as pd
from pathlib import Path
import tensorqtl
from tensorqtl import genotypeio, cis

parser = argparse.ArgumentParser(description='Run TensorQTL like the command line interface, but allow more options.')
parser.add_argument('--mode', help='Mapping mode')
parser.add_argument('geno_prefix', help='Genotype prefix')
parser.add_argument('bed', help='Input bed file')
parser.add_argument('out_prefix', help='Output prefix')
parser.add_argument('--covariates', help='Covariates file')
parser.add_argument('--cis_output', help='Output from cis mode')
parser.add_argument('--output_dir', help='Output directory')
parser.add_argument('--phenotype_groups', help='Header-less TSV with two columns: phenotype_id, group_id')
args = parser.parse_args()

# Copied from tensorqtl.py main():
if args.phenotype_groups is not None:
    group_s = pd.read_csv(args.phenotype_groups, sep='\t', index_col=0, header=None, squeeze=True)
    # verify sort order
    group_dict = group_s.to_dict()
    previous_group = ''
    parsed_groups = 0
    for i in phenotype_df.index:
        if group_dict[i]!=previous_group:
            parsed_groups += 1
            previous_group = group_dict[i]
    if not parsed_groups == len(group_s.unique()):
        raise ValueError('Groups defined in input do not match phenotype file (check sort order).')
else:
    group_s = None

# Actually I will just copy and modify the original tensorqtl.py.

if args.mode == 'cis':
    pheno, pheno_pos = tensorqtl.read_phenotype_bed(args.bed)
    covar = pd.read_csv(args.covariates, sep='\t', index_col=0).T
    pr = genotypeio.PlinkReader(args.geno_prefix)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    d = cis.map_cis(genotype_df, variant_df, pheno, pheno_pos, covar, random_tiebreak=True)
    tensorqtl.calculate_qvalues(d, qvalue_lambda=0.85)
    output = Path(args.output_dir) / f'{args.prefix}.cis_qtl.txt.gz'
    d.to_csv(output, sep='\t', float_format='%.6g')

elif args.mode == 'cis_independent':
    pheno, pheno_pos = tensorqtl.read_phenotype_bed(args.bed)
    covar = pd.read_csv(args.covariates, sep='\t', index_col=0).T
    pr = genotypeio.PlinkReader(args.geno_prefix)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    cis_df = pd.read_csv(args.cis_output, sep='\t', index_col=0)
    d = cis.map_independent(genotype_df, variant_df, cis_df, pheno, pheno_pos, covar, random_tiebreak=True)
    output = Path(args.output_dir) / f'{args.prefix}.cis_independent_qtl.txt.gz'
    d.to_csv(output, sep='\t', index=False, float_format='%.6g')

