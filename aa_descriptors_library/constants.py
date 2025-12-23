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

NON_PERMUTABLE_INDICES_SIDECHAIN = {
    "C": [2, 5, 6],
    "D": [2, 5, 6, 7, 8],
    "E": [2, 5, 8, 9, 10, 11],
    "F": [2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    "H": [2, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    "I": [2, 4, 7, 11],
    "K": [2, 5, 8, 11, 14],
    "L": [2, 5, 7, 11],
    "M": [2, 5, 8, 9],
    "N": [2, 5, 6, 7],
    "P": [2, 5, 8],
    "Q": [2, 5, 8, 9, 10],
    "R": [2, 5, 8, 11, 13, 14, 17],
    "S": [2, 5, 6],
    "T": [2, 4, 5, 6],
    "V": [2, 4, 8],
    "W": [2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
    "Y": [2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
}
