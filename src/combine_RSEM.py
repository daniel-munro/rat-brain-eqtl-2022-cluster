import argparse
import os
import gzip

parser = argparse.ArgumentParser(description="Combine gzipped RSEM output files into gzipped GCT file.")
parser.add_argument("infiles", nargs="+", help="List of RSEM files. Sample IDs will be part of filename after all '/' and before the next '.'.")
parser.add_argument("field", help="Which data column to extract.")
parser.add_argument("outfile", help="Name of output file (should end in '.gct.gz')")
args = parser.parse_args()

# infile = sys.argv[1]
# field = sys.argv[2]             # Which data column to use.
# RSEM columns: gene_id,transcript_id(s), length, effective_length,
# expected_count, TPM, FPKM

# with gzip.open(infile, "rt") as rsem:
#     lines = rsem.read().splitlines()
#     print("#1.2")
#     print("{}\t{}".format(len(lines) - 1, 1))
#     for i, line in enumerate(lines):
#         items = line.split("\t")
#         if i == 0:
#             # assert items[4] == "expected_count" and items[5] == "TPM"
#             datacol = items.index(field)
#         else:
#             # Just '.' for description column.
#             print("{}\t{}\t{}".format(items[0], ".", items[datacol]))

sample_ids = [os.path.split(x)[1].split('.')[0] for x in args.infiles]
rsem_tables = []              # List of tables, each is list of lists.
for infile in args.infiles:
    with gzip.open(infile, "rt") as rsem:
        lines = rsem.read().splitlines()
        items = [line.split("\t") for line in lines]
        rsem_tables.append(items)
data_col = rsem_tables[0][0].index(args.field)

with gzip.open(args.outfile, "wt", compresslevel=6) as f:
    f.write("#1.2\n")
    f.write("{}\t{}\n".format(len(rsem_tables[0]) - 1, len(rsem_tables)))
    f.write("{}\t{}\t{}\n".format("Name", "Description", "\t".join(sample_ids)))
    for i in range(1, len(rsem_tables[0])):
        gene = rsem_tables[0][i][0]
        values = "\t".join([table[i][data_col] for table in rsem_tables])
        f.write("{}\t{}\t{}\n".format(gene, ".", values))
