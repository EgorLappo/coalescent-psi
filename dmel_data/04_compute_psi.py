#!/usr/bin/env python 
import polars as pl
from tqdm import tqdm
import random
import os
import sys


def psi_quad(dd, i, eur_a, eur_b, afr_a, afr_b):
    dq = dd.select(
        [
            pl.col("boot"),
            pl.lit(i).alias("quad"),
            pl.col("fragment_id"),
            pl.col("site_id"),
            pl.col("SIM"),
            pl.col(eur_a).alias("eur_a"),
            pl.col(eur_b).alias("eur_b"),
            pl.col(afr_a).alias("afr_a"),
            pl.col(afr_b).alias("afr_b"),
        ]
    )

    return (
        dq.with_columns(
            pl.col("SIM")
            .list.concat("eur_a")
            .list.concat("eur_b")
            .list.concat("afr_a")
            .list.concat("afr_b")
            .alias("concat")
        )
        # filter deletions
        .filter(~pl.col("concat").list.contains("-"))
        # filter sites with >2 alleles
        .filter(pl.col("concat").list.n_unique() == 2)
        # mark derived allele
        .with_columns(
            pl.col("concat")
            .list.set_difference(pl.col("SIM").str.split(" "))
            .alias("derived")
            .list.get(0)
        )
        # count derived allele
        .with_columns(
            pl.col("eur_a")
            .list.concat("eur_b")
            .list.join("")
            .str.count_matches(pl.col("derived"))
            .cast(pl.Int32)
            .alias("psi_eur"),
            pl.col("afr_a")
            .list.concat("afr_b")
            .list.join("")
            .str.count_matches(pl.col("derived"))
            .cast(pl.Int32)
            .alias("psi_afr"),
        )
        # filter sites that are invariable in D. mel
        .filter(~((pl.col("psi_eur") == 2) & (pl.col("psi_afr") == 2)))
        # filter sites private to either D. mel population
        .filter(~((pl.col("psi_eur") == 0) | (pl.col("psi_afr") == 0)))
        # compute psi for remaining shared sites
        .with_columns((pl.col("psi_eur") - pl.col("psi_afr")).alias("psi"))
        .group_by(["boot", "quad", "fragment_id"])
        .agg(pl.col("psi"))
        .with_columns(
            pl.col("psi").list.sample(1).list.get(0),
            pl.lit(",".join(sorted([eur_a, eur_b]))).alias("eur_lines"),
            pl.lit(",".join(sorted([afr_a, afr_b]))).alias("afr_lines"),
        )
        .collect()
        .sample(1)
    )


d = pl.scan_parquet("data/sites_boot.parquet")

lines = [c for c in d.collect_schema().names() if "MEL" in c]
afr_lines = [ln for ln in lines if "ZBMEL" in ln]
eur_lines = [ln for ln in lines if "MEL" in ln and "ZB" not in ln]

iboot = int(sys.argv[1])
print(f"working on bootstrap replicate {iboot}")
dd = d.filter(pl.col("boot") == iboot)

dboots = []
quads = []

for i in tqdm(range(1000), total=1000):
    eur_a, eur_b = random.sample(eur_lines, k=2)
    afr_a, afr_b = random.sample(afr_lines, k=2)
    dboots.append(psi_quad(dd, i, eur_a, eur_b, afr_a, afr_b))

dboot = pl.concat(dboots)

os.makedirs("data/boots", exist_ok=True)
dboot.write_parquet(f"data/boots/boot_{iboot}.parquet")
