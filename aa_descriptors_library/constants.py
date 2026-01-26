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

BACKBONE_SMARTS = {
    "zwitterion": "[NH3+][CH]({R})[C](=[O])[O-]",
    "neutral": "[NH2][CH]({R})[C](=[O])[OH]",
    "pro_zwitterion": "[NH2+]1{R}[CH]1[C](=[O])[O-]",
    "pro_neutral": "[NH]1{R}[CH]1[C](=[O])[OH]",
}

SIDECHAIN_SMARTS = {
    "C": {0: "[CH2][SH]", -1: "[CH2][S-]"},
    "D": {0: "[CH2][C](=[O])[OH]", -1: "[CH2][C](=[O])[O-]"},
    "E": {0: "[CH2][CH2][C](=[O])[OH]", -1: "[CH2][CH2][C](=[O])[O-]"},
    "F": {0: "[CH2][c]1[cH][cH][cH][cH][cH]1"},
    "H": {
        0: {"D": "[CH2][c]1[nH][cH][n][cH]1", "E": "[CH2][c]1[n][cH][nH][cH]1"},
        -1: "[CH2][c]1[n][cH][n][cH]1",
        +1: "[CH2][c]1[nH][cH][nH][cH]1",
    },
    "I": {0: "[CH]([CH3])[CH2][CH3]"},
    "K": {0: "[CH2][CH2][CH2][CH2][NH2]", +1: "[CH2][CH2][CH2][CH2][NH3+]"},
    "L": {0: "[CH2][CH]([CH3])[CH3]"},
    "M": {0: "[CH2][CH2][S][CH3]"},
    "N": {0: "[CH2][C](=[O])[NH2]"},
    "P": {0: "[CH2][CH2][CH2]"},
    "Q": {0: "[CH2][CH2][C](=[O])[NH2]"},
    "R": {+1: "[CH2][CH2][CH2][NH][C](=[NH2+])[NH2]"},
    "S": {0: "[CH2][OH]"},
    "T": {0: "[CH]([CH3])[OH]"},
    "V": {0: "[CH]([CH3])[CH3]"},
    "W": {0: "[CH2][c]1[cH][nH][c]2[cH][cH][cH][cH][c]12"},
    "Y": {0: "[CH2][c]1[cH][cH][c]([OH])[cH][cH]1"},
}

OXT_IDX = {
    "C": {0: 14, -1: 13},
    "D": {0: 16, -1: 15},
    "E": {0: 19, -1: 18},
    "F": {0: 23},
    "H": {0: {"D": 20, "E": 20}, -1: 19, +1: 21},
    "I": {0: 22},
    "K": {0: 24, +1: 25},
    "L": {0: 22},
    "M": {0: 20},
    "N": {0: 17},
    "P": {0: 17},
    "Q": {0: 20},
    "R": {+1: 27},
    "S": {0: 14},
    "T": {0: 17},
    "V": {0: 19},
    "W": {0: 27},
    "Y": {0: 24},
}
