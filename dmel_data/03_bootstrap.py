#!/usr/bin/env python
import os
import random
from tqdm import tqdm
import pandas as pd
from Bio import AlignIO
from collections import defaultdict

values = defaultdict(dict)

fids = set()

for fragment in [
    f"data/fragments_aligned/{f}" for f in os.listdir("data/fragments_aligned")
]:
    seqs = AlignIO.read(fragment, "fasta")

    fid = int(fragment.split("/")[-1].split(".")[0])
    fids.add(fid)

    for s in seqs:
        line = s.id.split("_")[-1]
        values[line][fid] = s

afr_lines = [ln for ln in values.keys() if "ZBMEL" in ln]
eur_lines = [ln for ln in values.keys() if "MEL" in ln and "ZB" not in ln]

dfgs = []

sim_fragments = values["SIM"]
for fid in fids:
    rows = []

    for i in range(len(sim_fragments[fid])):
        # row = {k: values[k][fid][i] for k in values.keys()}
        row = {}
        for k in values.keys():
            if fid in values[k]:
                row[k] = values[k][fid][i]
        row["fragment_id"] = fid
        row["site_id"] = i
        rows.append(row)

    dfg = pd.DataFrame.from_records(rows)
    dfgs.append(dfg)

dg = pd.concat(dfgs)

dg.to_parquet("data/sites.parquet")

boots = []

fragments = list(dg["fragment_id"].unique())
for i in tqdm(range(1000), total=1000):
    boot_fragments = random.sample(fragments, len(fragments))

    for frag in boot_fragments:
        dfrag = dg.query("fragment_id == @frag").copy()
        dfrag["boot"] = i
        boots.append(dfrag)

boots = pd.concat(boots)

boots.to_parquet("data/sites_boot.parquet", index=False)
