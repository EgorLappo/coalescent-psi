#!/usr/bin/env python
import os
import subprocess
import pandas as pd

MUSCLE_EXE = "/Users/egor/miniforge3/envs/muscle/bin/muscle"

# first, make fasta files
os.makedirs("data/fragments", exist_ok=True)

# then, write
d = pd.read_csv("data/sequences.csv").reset_index().rename(columns={"index": "seqid"})

fragments = d["fragment_id"].unique()
for fid in fragments:
    dfrag = d[d["fragment_id"] == fid]

    with open(f"data/fragments/{fid}.fasta", "w") as f:
        for i, row in dfrag.iterrows():
            f.write(f">{row['sample']}_{row['fragment_id']}_{row['line_id']}\n")
            f.write(f"{row['seq']}\n")

# now run muscle
os.makedirs("data/fragments_aligned", exist_ok=True)

for fid in fragments:
    cmd = f"{MUSCLE_EXE} -align data/fragments/{fid}.fasta -output data/fragments_aligned/{fid}.fasta"
    subprocess.run(cmd, shell=True, check=True)