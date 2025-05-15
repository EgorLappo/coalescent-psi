#!/usr/bin/env python
import os
import pandas as pd
from Bio import SeqIO


""" AJ ACCESSIONS """

AJ_accessions = ["data/embl/" + f for f in os.listdir("data/embl") if "AJ" in f]
records = [SeqIO.read(f, "embl") for f in AJ_accessions]

d = pd.DataFrame(
    {
        "id": [r.id for r in records],
        "seq": [str(r.seq) for r in records],
        "seqlen": [len(r.seq) for r in records],
        "feature_types": [" ".join([f.type for f in r.features]) for r in records],
    }
)
d["kind"] = d["feature_types"].apply(
    lambda x: "intron" if "intron" in x else "intergenic"
)
d["fragment"] = ""
d["sample"] = ""
d["chrom"] = ""
d["organism"] = ""
d = d.drop(columns=["feature_types"])

for i, row in d.iterrows():
    # deal with ids
    if row["kind"] == "intron":
        f = [feat for feat in records[i].features if feat.type == "intron"][0]
        d.at[i, "fragment"] = f.qualifiers["note"][0].strip("fragment ID ")
    elif row["kind"] == "intergenic":
        f = [feat for feat in records[i].features if feat.type == "misc_feature"][0]
        d.at[i, "fragment"] = f.qualifiers["note"][0].strip("fragment ID ")

    # other info
    f = [feat for feat in records[i].features if feat.type == "source"][0]
    if "country" in f.qualifiers:
        d.at[i, "sample"] = f.qualifiers["country"][0].lower()
    else:
        d.at[i, "sample"] = "simulans"
    d.at[i, "chrom"] = f.qualifiers["chromosome"][0]
    d.at[i, "organism"] = f.qualifiers["organism"][0]


d["fragment"] = d["fragment"].str.upper()
dmel = d[d["organism"] == "Drosophila melanogaster"].copy()
dsim = d[d["organism"] == "Drosophila simulans"].copy()

# for some reason, there is an id like 272sim for two D. mel samples.... (AJ570016 and AJ570044)
# also african populations are labeled as ZBMEL...
dmel["fragment_id"] = dmel["fragment"].apply(
    lambda x: x.split("SIM")[0]
    if ("SIM" in x)
    else x.split("MEL")[0].strip("_").strip("ZB")
)

# here we also fill in for the two samples with 'SIM' in the fragment name
dmel["line_id"] = dmel["fragment"].apply(
    lambda x: "MEL01" if "SIM" in x else x.lstrip("0123456789.- ")
)  # .apply(lambda x: '01' if len(x.split('MEL')) == 1 else x.split('MEL')[1].strip("_"))

# parallel thing happens with D. simulans but in the opposite direction
dsim["fragment_id"] = dsim["fragment"].apply(
    lambda x: x.split("MEL")[0] if ("MEL" in x) else x.split("SIM")[0].strip("_")
)

d_aj = pd.concat([dmel, dsim])
d_aj["fragment_id"] = d_aj["fragment_id"].apply(lambda x: int(x))
d_aj["line_id"] = d_aj["line_id"].fillna("SIM")

""" AM ACCESSIONS """

AM_accessions = ["data/embl/" + f for f in os.listdir("data/embl") if "AM" in f]
records = [SeqIO.read(f, "embl") for f in AM_accessions]

d = pd.DataFrame(
    {
        "id": [r.id for r in records],
        "seq": [str(r.seq) for r in records],
        "seqlen": [len(r.seq) for r in records],
        "feature_types": [" ".join([f.type for f in r.features]) for r in records],
    }
)
d["kind"] = d["feature_types"].apply(lambda x: "intron" if "intron" in x else "STS")
d["fragment_id"] = ""
d["sample"] = ""
d["chrom"] = ""
d["organism"] = ""
d["line_id"] = ""
d = d.drop(columns=["feature_types"])

for i, row in d.iterrows():
    # get info from features
    if row["kind"] == "intron":
        f = [feat for feat in records[i].features if feat.type == "intron"][0]
        f_source = [feat for feat in records[i].features if feat.type == "source"][0]

        d.at[i, "fragment_id"] = f_source.qualifiers["clone"][0].strip("fragment ID ")
    elif row["kind"] == "STS":
        f = [feat for feat in records[i].features if feat.type == "STS"][0]
        d.at[i, "fragment_id"] = f.qualifiers["standard_name"][0].strip("fragment ID ")

    # other info
    f = [feat for feat in records[i].features if feat.type == "source"][0]
    if "note" in f.qualifiers:
        d.at[i, "line_id"] = f.qualifiers["note"][0].strip("line: ").upper()

    if "country" in f.qualifiers:
        d.at[i, "sample"] = f.qualifiers["country"][0].lower()
    else:
        d.at[i, "sample"] = "simulans"
    d.at[i, "chrom"] = f.qualifiers["chromosome"][0]
    d.at[i, "organism"] = f.qualifiers["organism"][0]

# replace sample names
d["sample"] = d["sample"].replace(
    {"netherlands:leiden": "netherlands", "usa:ca, davis": "simulans"}
)
d_am = d.copy()

# fix line id issues
d_am["line_id"] = d_am["line_id"].apply(lambda x: "SIM" if x == "" else x)
d_am["line_id"] = d_am["line_id"].apply(lambda x: "MEL01" if x == "MEL1" else x)
d_am["line_id"] = d_am["line_id"].apply(lambda x: "MEL02" if x == "MEL2" else x)
d_am["line_id"] = d_am["line_id"].apply(lambda x: "MEL20" if x == "MELL20" else x)
d_am["line_id"] = d_am["line_id"].apply(lambda x: "MEL13" if x == "MEL13_ALIGN" else x)
d_am["line_id"] = d_am["line_id"].apply(lambda x: "ZBMEL82" if x == "ZBMEL82_N" else x)

""" COMBINE """

d = pd.concat([d_aj, d_am])

seq_counts = (
    d.groupby(["fragment_id", "sample"])
    .agg({"seq": "count"})
    .reset_index()
    .pivot(index="fragment_id", columns="sample", values="seq")
)

seq_counts.to_csv("data/seq_counts.csv")

# drop fragments for which we don't have data in all three populations
d = d[d["fragment_id"].isin(seq_counts.dropna().index)]

d.to_csv("data/sequences.csv", index=False)