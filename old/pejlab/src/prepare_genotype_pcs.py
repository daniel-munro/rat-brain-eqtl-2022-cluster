import sys

evec = sys.argv[1]
outfile = sys.argv[2]
N_PC = int(sys.argv[3])

lines = open(evec, "r").read().splitlines()[1:]  # Eigenvalues are in first line.
lines = [line.split() for line in lines]
rows = list(zip(*lines))[:(N_PC + 1)]
rows[0] = [ind.replace("0:", "") for ind in rows[0]]
index = ["ID"] + ["PC" + str(i + 1) for i in range(N_PC)]

with open(outfile, "w") as out:
    for i, row in enumerate(rows):
        out.write("{}\t{}\n".format(index[i], "\t".join(row)))
