import argparse
from gtfparse import read_gtf
import numpy as np
from pathlib import Path
import pandas as pd


def load_tss(anno: pd.DataFrame) -> pd.DataFrame:
    anno = anno.loc[anno['feature'] == 'gene', :]
    anno['chromEnd'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['chromStart'] = anno['chromEnd'] - 1  # BED coordinates are 0-based
    # anno = anno.rename(columns={'seqname': '#chrom', 'gene_id': 'name'})
    anno['#chrom'] = anno['seqname']
    anno = anno.sort_values(['#chrom', 'chromStart'])
    return anno[['#chrom', 'chromStart', 'chromEnd', 'gene_id']]


def load_exons(anno: pd.DataFrame) -> pd.DataFrame:
    anno = anno.loc[anno['feature'] == 'exon', :]
    anno['chrom'] = anno['seqname']
    anno['exonStart'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['exonEnd'] = np.where(anno['strand'] == '+', anno['end'], anno['start'])
    return anno[['gene_id', 'chrom', 'exonStart', 'exonEnd']]


def map_introns_to_genes(introns: list, exons: pd.DataFrame) -> pd.DataFrame:
    exons['exonStart'] = exons['exonStart'].astype(str)
    exons['exonEnd'] = exons['exonEnd'].astype(str)
    df = pd.DataFrame({'intron': introns})
    df[['chrom', 'chr_start', 'chr_end', 'clu', 'cluster', 'strand']] = df['intron'].str.split(r':|_', expand=True)
    df['start'] = np.where(df['strand'] == '+', df['chr_start'], df['chr_end'])
    df['end'] = np.where(df['strand'] == '+', df['chr_end'], df['chr_start'])
    start_matches = df.merge(
        exons[['chrom', 'exonEnd', 'gene_id']].rename(columns={'exonEnd': 'start'}),
        on=['chrom', 'start'],
        how='inner'
    )
    end_matches = df.merge(
        exons[['chrom', 'exonStart', 'gene_id']].rename(columns={'exonStart': 'end'}),
        on=['chrom', 'end'],
        how='inner'
    )
    df = pd.concat([start_matches, end_matches])
    clust_genes = df.groupby('cluster', group_keys=False).agg({'gene_id': pd.Series.unique})
    clust_genes = clust_genes.reset_index().explode('gene_id')
    return clust_genes


parser = argparse.ArgumentParser(description='Create BED file with splicing phenotypes')
parser.add_argument('input', type=Path, help='LeafCutter cluster file (*_perind_numers.counts.gz)')
parser.add_argument('gtf', type=Path, help="GTF annotation file")
parser.add_argument('output', type=Path, help="Output path (*.bed)")
args = parser.parse_args()

df = pd.read_csv(args.input, sep=' ')
samples = list(df.columns)
df.index = df.index.rename('intron')
df = df.reset_index()
df['cluster'] = df['intron'].str.extract(r'clu_(\d+)_', expand=False)
for col in samples:
    df[col] = df.groupby('cluster', group_keys=False).apply(lambda g: g[col] / g[col].sum())
anno = read_gtf(args.gtf)
exons = load_exons(anno)
genes = map_introns_to_genes(df['intron'], exons)
df = df.merge(genes, on='cluster', how='left')
tss = load_tss(anno)
df = tss.merge(df, on='gene_id', how='inner')
df['name'] = df['gene_id'] + ':' + df['intron']
df = df[['#chrom', 'chromStart', 'chromEnd', 'name'] + samples]
df.to_csv(args.output, sep='\t', index=False, float_format='%g')
