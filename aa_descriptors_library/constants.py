SQL_PATH = "/cluster/project/jorner/lajacot/projects/aa-descriptors-library/data/merged_whole_filtered_library_sql.db"

NDRD_PATH = (
    "/cluster/project/jorner/lajacot/projects/aa-descriptors-library/data/ndrd.csv"
)

NUMBER_OF_CHI_ANGLES = {
    "C": 1,
    "S": 1,
    "T": 1,
    "V": 1,
    "N": 2,
    "D": 2,
    "H": 2,
    "I": 2,
    "L": 2,
    "F": 2,
    "P": 2,
    "W": 2,
    "Y": 2,
    "Q": 3,
    "E": 3,
    "M": 3,
    "R": 4,
    "K": 4,
}

ONE_TO_THREE_AA = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

THREE_TO_ONE_AA = {v: k for k, v in ONE_TO_THREE_AA.items()}
