import sys
import pysam

INFILE = sys.argv[1]
OUTFILE = sys.argv[2]

vcf = pysam.VariantFile(INFILE)

with open(OUTFILE, "w") as out:
    out.write("marker,chr,pos\n")
    for rec in vcf.fetch():
        ID = rec.id if rec.id is not None else "{}:{}".format(rec.contig, rec.pos)
        pos = rec.pos / 1000000  # Units are Mbp.
        out.write("{},{},{}\n".format(ID, rec.contig, pos))
