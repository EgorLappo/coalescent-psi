#!/usr/bin/env python

import os
import requests
from tqdm import tqdm

# ids:
#     * EMBL: AJ568984 to AJ571588,  AM000058 to AM003900



if __name__ == "__main__":
    os.makedirs("data/embl", exist_ok=True)

    EMBL = [f"AJ{i}" for i in range(568984, 571589)]
    EMBL += [f"AM{i:06d}" for i in range(58, 3901)]

    # dont download files that already exist
    EMBL = [eid for eid in EMBL if not os.path.exists(f"data/embl/{eid}.embl")]

    for eid in tqdm(EMBL, desc="Downloading EMBL records"):
        # url of the form https://www.ebi.ac.uk/ena/browser/api/embl/AJ568984
        url = f"https://www.ebi.ac.uk/ena/browser/api/embl/{eid}"
        response = requests.get(url)
        filename = f"data/embl/{eid}.embl"
        with open(f"data/embl/{eid}.embl", "w") as f:
            f.write(response.text)

