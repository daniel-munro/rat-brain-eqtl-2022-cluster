import sys
import numpy as np
import pandas as pd
# from rpy2 import robjects
# from rpy2.robjects import numpy2ri

WINDOW = 50

GENO_FILE = sys.argv[1]
FOUNDER_FILE = sys.argv[2]
OUTFILE = sys.argv[3]
N_SNPS = int(sys.argv[4])

# geno = pd.read_csv("data/qtl2/chr/geno_chr12.csv", index_col="id").to_numpy()
# founder = pd.read_csv("data/qtl2/chr/founder_geno_chr12.csv", index_col="id").to_numpy()
geno = pd.read_csv(GENO_FILE, index_col="id").to_numpy()
founder = pd.read_csv(FOUNDER_FILE, index_col="id").to_numpy()

founder[founder == "A"] = 0
founder[founder == "B"] = 1
founder[founder == "-"] = -10
founder = founder.astype(int)

geno[geno == "A"] = 0
geno[geno == "H"] = 1
geno[geno == "B"] = 2
geno[geno == "-"] = -100
geno = geno.astype(int)

hap_pairs = np.zeros((founder.shape[0], 36), dtype=int)

pair_indices = []
col = 0
for i in range(7):
    for j in range(i + 1, 8):
        hap_pairs[:, col] = founder[:, i] + founder[:, j]
        pair_indices.append((i, j))
        col += 1
for i in range(8):
    hap_pairs[:, col] = founder[:, i] + founder[:, i]
    pair_indices.append((i, j))
    col += 1

similarity = np.zeros((geno.shape[1], hap_pairs.shape[1], geno.shape[0]))
# similarity = np.zeros((10, hap_pairs.shape[1], geno.shape[0]))
for sample in range(geno.shape[1]):
    for hap_pair in range(hap_pairs.shape[1]):
        for snp in range(geno.shape[0]):
            snp_left = max(0, snp - WINDOW)
            snp_right = min(snp + WINDOW + 1, geno.shape[0] - 1)
            genos1 = geno[snp_left:(snp_right), sample]
            genos2 = hap_pairs[snp_left:(snp_right), hap_pair]
            similarity[sample, hap_pair, snp] = np.mean(genos1 == genos2)

# collapsed = np.zeros((similarity.shape[0], 8, similarity.shape[2]))
# for i in range(similarity.shape[1]):
#     found1, found2 = pair_indices[i]
#     collapsed[:, found1, :] += similarity[:, i, :] / 2
#     collapsed[:, found2, :] += similarity[:, i, :] / 2

top_pair = np.argmax(similarity, axis=1)
collapsed = np.zeros((top_pair.shape[0], 8, top_pair.shape[1]), dtype=int)
for i in range(top_pair.shape[0]):
    for j in range(top_pair.shape[1]):
        found1, found2 = pair_indices[top_pair[i, j]]
        collapsed[i, found1, j] += 1
        collapsed[i, found2, j] += 1

# Save probs for only every Nth  (e.g. to get 2000 loci per chromosome) for plotting, etc.
subset = [int(i * (collapsed.shape[2] / 2000)) for i in range(2000)]
collapsed = collapsed[:, :, subset]

# saveRDS = robjects.r["saveRDS"]
# saveRDS(similarity, "analysis/haplotype_sim_chr12.rds")
# np.save("analysis/haplotype_sim_chr12.npy", collapsed)
np.save(OUTFILE, collapsed)
