"""
Create the single-rotamer entries for alanine and glycine (which have no chi angles
and are therefore absent from the Dunbrack rotamer library), and combine them with
the Dunbrack library data into the full angles file used to generate the database.
"""

import pandas as pd

from rotadel.common.config import ROOT_DIR

single_aas_json = ROOT_DIR / "data" / "angles_single_AAs.json"
dunbrack_json = ROOT_DIR / "data" / "angles_Dunbrack.json"
combined_json = ROOT_DIR / "data" / "angles_combined.json"

# Alanine and glycine have no side-chain dihedral angles => single "rotamer" each, with probability 1.0
single_aas_data = {
    "Aa0a0r0000": {
        "res": "ALA",
        "letter": "A",
        "phi": None,
        "psi": None,
        "chi1": None,
        "chi2": None,
        "chi3": None,
        "chi4": None,
        "prob": 1.0,
    },
    "Ga0a0r0000": {
        "res": "GLY",
        "letter": "G",
        "phi": None,
        "psi": None,
        "chi1": None,
        "chi2": None,
        "chi3": None,
        "chi4": None,
        "prob": 1.0,
    },
}
single_aas_df = pd.DataFrame.from_dict(single_aas_data, orient="index")
single_aas_df.to_json(single_aas_json, orient="index", indent=4)

# Combine into angles_combined.json
with open(single_aas_json, "rb") as f:
    single_aas_bytes = f.read()
with open(dunbrack_json, "rb") as f:
    dunbrack_bytes = f.read()
single_aas_body = single_aas_bytes[:-2]  # strip trailing "\n}"
dunbrack_body = dunbrack_bytes[dunbrack_bytes.index(b"{") + 1 :]  # noqa: E203
combined_bytes = single_aas_body + b"," + dunbrack_body
with open(combined_json, "wb") as f:
    f.write(combined_bytes)
