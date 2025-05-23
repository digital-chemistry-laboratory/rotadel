"""
Save the rotamers data from Dunbrack library into a JSON file.
The angles and probabilities are extracted, and a unique ID is generated for each rotamer.
"""

import pandas as pd

lib_file = "/cluster/project/jorner/lajacot/projects/aa-descriptors-library/Dunbrack_library/ALL_rotamers.lib"
json_output = "/cluster/project/jorner/lajacot/projects/aa-descriptors-library/gen_lib/data/angles_Dunbrack.json"

column_names = [
    "res",
    "phi",
    "psi",
    "count",
    "r1",
    "r2",
    "r3",
    "r4",
    "prob",
    "chi1",
    "chi2",
    "chi3",
    "chi4",
    "chi1Sig",
    "chi2Sig",
    "chi3Sig",
    "chi4Sig",
]
lib_df = pd.read_csv(
    lib_file,
    sep=r"\s+",
    comment="#",
    names=column_names,
    low_memory=False,
)

# Drop the TPR, CPR, CYH and CYD rows
# TPR/CPR are trans/cis prolines and both already included in PRO
# CYH/CYD are nondisulfide/disulfide-bonded cysteines and both already included in CYS
lib_df = lib_df[~lib_df["res"].isin(["TPR", "CPR", "CYH", "CYD"])]

amino_acid_codes = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
lib_df["letter"] = lib_df["res"].map(amino_acid_codes)

# Generate rotamer ID
lib_df["rotamer_id"] = (
    lib_df["letter"]
    + ":"
    + lib_df["phi"].astype(str)
    + ":"
    + lib_df["psi"].astype(str)
    + ":"
    + "r"
    + lib_df[["r1", "r2", "r3", "r4"]].astype(str).agg("".join, axis=1)
)
lib_df.set_index("rotamer_id", inplace=True)

lib_df.drop(
    ["count", "r1", "r2", "r3", "r4", "chi1Sig", "chi2Sig", "chi3Sig", "chi4Sig"],
    axis="columns",
    inplace=True,
)
ordered_columns_names = [
    "res",
    "letter",
    "phi",
    "psi",
    "chi1",
    "chi2",
    "chi3",
    "chi4",
    "prob",
]
lib_df = lib_df[ordered_columns_names]

# The chi angles with 0.0 should be NaN. Replacing them by None in pandas will then write null in json.
for col in ["chi1", "chi2", "chi3", "chi4"]:
    lib_df[col] = lib_df[col].where(lib_df[col] != 0.0, None)

lib_df.to_json(
    json_output,
    orient="index",
    indent=4,
)
