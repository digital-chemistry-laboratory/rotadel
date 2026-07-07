import json
import pandas as pd

from rotadel.common.config import ROOT_DIR

"""
Create the rotamer IDs containing the charge and tautomer information.

The ID format is (all attached):
    X: one letter code of the amino acid
    axxxx: phi angle
        'a': indicates an angle
        then max 4 symbols for the angle value: '-' sign if negative and 2-3 digits
    axxxx: psi angle
        same as above
    rxxxx: numerical designation of the rotamer
        'r': indicates rotamer designation from Dunbrack's library
        then 4 numbers for r1, r2, r3, r4, as reported in Dunbrack's library
    cxx: charge of the species
        'c': indicates a charge
        then either '-1', '0', '+1' for the charge
    tX: tautomeric species, only for histidine
        't': indicates tautomeric information
        then either 'D' or 'E' for which histidine nitrogen is protonated
    Example: Ha-180a40r11200c0tE

Save into a csv file to use to with GNU parallel to run calculations on all rotamers.
"""

dunbrack_json = ROOT_DIR / "data" / "angles_Dunbrack.json"
with open(dunbrack_json, "r") as f:
    all_rotamers_data = json.load(f)

ids_species = []
for rotamer_id, entry in all_rotamers_data.items():
    letter = entry["letter"]
    if letter in ["D", "E", "C"]:
        ids_species.append(f"{rotamer_id}c-1")
        ids_species.append(f"{rotamer_id}c0")
    elif letter in ["H"]:
        ids_species.append(f"{rotamer_id}c0tD")
        ids_species.append(f"{rotamer_id}c0tE")
        ids_species.append(f"{rotamer_id}c+1")
        ids_species.append(f"{rotamer_id}c-1")
    elif letter in ["R"]:
        ids_species.append(f"{rotamer_id}c+1")
    elif letter in ["K"]:
        ids_species.append(f"{rotamer_id}c+1")
        ids_species.append(f"{rotamer_id}c0")
    else:
        ids_species.append(f"{rotamer_id}c0")
df = pd.DataFrame(ids_species)

df.to_csv("data/keys_species_initial.csv", index=False, header=False)
