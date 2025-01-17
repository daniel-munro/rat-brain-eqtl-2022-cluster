import pandas as pd

# lm_files = [f"data/gemma/lm/Acbc.{chrn}.assoc.txt" for
# modifier = "_cov"
modifier = ""

pairs = []
for chrn in range(1, 21):
    print(chrn)
    lm_file = f"data/gemma/lm{modifier}/Acbc.{chrn}.assoc.txt"
    # lmm_file = f"data/gemma/lmm{modifier}/Acbc.{chrn}.assoc.txt"
    lmm_file = f"data/gemma/lm_cov/Acbc.{chrn}.assoc.txt"
    # round_trip prevents precision issues in output
    d = pd.read_csv(lm_file, sep="\t", float_precision="round_trip")
    d = d.groupby("gene_id").sample(1)
    d = d[["gene_id", "rs", "p_wald"]]
    d = d.rename(columns={"p_wald": "p_lm"})

    d2 = pd.read_csv(lmm_file, sep="\t", float_precision="round_trip")
    d2 = d2[["gene_id", "rs", "p_wald"]]
    d2 = d2.rename(columns={"p_wald": "p_lm"})

    d = d.merge(d2, on=["gene_id", "rs"], how="left")
    pairs.append(d)

pd.concat(pairs).to_csv(f"data/gemma/Acbc{modifier}_cov2.random_pairs.txt",
                        sep="\t", index=False)
