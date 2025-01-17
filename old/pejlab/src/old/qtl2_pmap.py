import sys
# import pysam

INFILE = sys.argv[1]
OUTFILE = sys.argv[2]

# vcf = pysam.VariantFile(INFILE)
IDs = open(INFILE, "r").read().splitlines()

with open(OUTFILE, "w") as out:
    out.write("marker,chr,pos\n")
    # for rec in vcf.fetch():
    for ID in IDs:
        # ID = rec.id if rec.id is not None else "{}:{}".format(rec.contig, rec.pos)
        chrom, pos = tuple(ID.replace("chr", "").split(":"))
        # pos = rec.pos / 1000000  # Units are Mbp.
        pos = int(pos) / 1e6  # Units are Mbp.
        # out.write("{},{},{}\n".format(ID, rec.contig, pos))
        out.write("{},{},{}\n".format(ID, chrom, pos))
