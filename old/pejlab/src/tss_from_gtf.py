import argparse
import pandas as pd
import numpy as np
from gtfparse import read_gtf

parser = argparse.ArgumentParser(description="Extract gene transcription start sites from GTF file.")
parser.add_argument("infile", help="GTF annotation file.")
parser.add_argument("outfile", help="Output file name ('*.tsv').")
args = parser.parse_args()

# anno = read_gtf("data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf")
anno = read_gtf(args.infile)
anno = anno[anno["feature"] == "gene"]
anno["start"][anno["strand"] == "-"] = anno["end"][anno["strand"] == "-"]
anno = anno[anno["seqname"].isin([str(i) for i in range(1, 21)])]
anno["seqname"] = anno["seqname"].astype(int)  # for sorting
anno = anno.rename(columns={"seqname": "chr", "start": "pos"})
anno = anno[["gene_id", "chr", "pos"]]
anno = anno.sort_values(["chr", "pos"])

anno.to_csv(args.outfile, sep="\t", index=False)#, float_format="%g")
