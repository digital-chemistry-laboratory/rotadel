NUMBER_OF_CHI_ANGLES = {
    "A": 0,
    "G": 0,
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
    "A": ["CB"],
    "C": ["CB", "SG", "HG"],
    "D": ["CB", "CG", "OD1", "OD2", "HD2"],
    "E": ["CB", "CG", "CD", "OE1", "OE2", "HE2"],
    "F": [
        "CB",
        "CG",
        "CD1",
        "HD1",
        "CE1",
        "HE1",
        "CD2",
        "HD2",
        "CE2",
        "HE2",
        "CZ",
        "HZ",
    ],
    "G": [],
    "H": ["CB", "CG", "ND1", "HD1", "CD2", "HD2", "CE1", "HE1", "NE2", "HE2"],
    "I": ["CB", "HB", "CG1", "CG2", "CD1"],
    "K": ["CB", "CG", "CD", "CE", "NZ"],
    "L": ["CB", "CG", "HG", "CD1", "CD2"],
    "M": ["CB", "CG", "SD", "CE"],
    "N": ["CB", "CG", "OD1", "ND2"],
    "P": ["CB", "CG", "CD"],
    "Q": ["CB", "CG", "CD", "OE1", "NE2"],
    "R": ["CB", "CG", "CD", "NE", "HE", "CZ", "NH1", "NH2"],
    "S": ["CB", "OG", "HG"],
    "T": ["CB", "HB", "OG1", "HG1", "CG2"],
    "V": ["CB", "HB", "CG1", "CG2"],
    "W": [
        "CB",
        "CG",
        "CD1",
        "HD1",
        "CD2",
        "NE1",
        "HE1",
        "CE2",
        "CE3",
        "HE3",
        "CZ2",
        "HZ2",
        "CZ3",
        "HZ3",
        "CH2",
        "HH2",
    ],
    "Y": [
        "CB",
        "CG",
        "CD1",
        "HD1",
        "CE1",
        "HE1",
        "CD2",
        "HD2",
        "CE2",
        "HE2",
        "CZ",
        "OH",
        "HH",
    ],
}

BACKBONE_SMARTS = {
    "zwitterion": "[NH3+][CH]({R})[C](=[O])[O-]",
    "neutral": "[NH2][CH]({R})[C](=[O])[OH]",
    "pro_zwitterion": "[NH2+]1{R}[CH]1[C](=[O])[O-]",
    "pro_neutral": "[NH]1{R}[CH]1[C](=[O])[OH]",
    "gly_zwitterion": "[NH3+][CH2][C](=[O])[O-]",
    "gly_neutral": "[NH2][CH2][C](=[O])[OH]",
}

SIDECHAIN_SMARTS = {
    "A": {0: "[CH3]"},
    "C": {0: "[CH2][SH]", -1: "[CH2][S-]"},
    "D": {0: "[CH2][C](=[O])[OH]", -1: "[CH2][C](=[O])[O-]"},
    "E": {0: "[CH2][CH2][C](=[O])[OH]", -1: "[CH2][CH2][C](=[O])[O-]"},
    "F": {0: "[CH2][c]1[cH][cH][cH][cH][cH]1"},
    "G": {0: "[H]"},
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

SIDECHAIN_PKA = {
    "D": 3.65,
    "E": 4.25,
    "H": 6.0,
    "C": 8.18,
    "K": 10.53,
}

RES_SPECIES = {
    # Alanine
    "ALA": {"letter": "A", "charge": 0, "tautomer": None},
    # Arginine
    "ARG": {"letter": "R", "charge": +1, "tautomer": None},
    # Asparagine
    "ASN": {"letter": "N", "charge": 0, "tautomer": None},
    # Aspartic acid
    "ASP": {"letter": "D", "charge": -1, "tautomer": None},
    "ASH": {"letter": "D", "charge": 0, "tautomer": None},
    "ASPP": {"letter": "D", "charge": 0, "tautomer": None},
    # Cysteine
    "CYS": {"letter": "C", "charge": 0, "tautomer": None},
    "CYX": {"letter": "C", "charge": 0, "tautomer": None},
    "CYM": {"letter": "C", "charge": -1, "tautomer": None},
    # Glutamic acid
    "GLU": {"letter": "E", "charge": -1, "tautomer": None},
    "GLH": {"letter": "E", "charge": 0, "tautomer": None},
    "GLUP": {"letter": "E", "charge": 0, "tautomer": None},
    # Glutamine
    "GLN": {"letter": "Q", "charge": 0, "tautomer": None},
    # Glycine
    "GLY": {"letter": "G", "charge": 0, "tautomer": None},
    # Histidine
    "HIS": {"letter": "H", "charge": 0, "tautomer": "E"},
    "HIE": {"letter": "H", "charge": 0, "tautomer": "E"},
    "HSE": {"letter": "H", "charge": 0, "tautomer": "E"},
    "HID": {"letter": "H", "charge": 0, "tautomer": "D"},
    "HSD": {"letter": "H", "charge": 0, "tautomer": "D"},
    "HIP": {"letter": "H", "charge": +1, "tautomer": None},
    "HIH": {"letter": "H", "charge": +1, "tautomer": None},
    "HSP": {"letter": "H", "charge": +1, "tautomer": None},
    # Isoleucine
    "ILE": {"letter": "I", "charge": 0, "tautomer": None},
    # Leucine
    "LEU": {"letter": "L", "charge": 0, "tautomer": None},
    # Lysine
    "LYS": {"letter": "K", "charge": +1, "tautomer": None},
    "LSN": {"letter": "K", "charge": 0, "tautomer": None},
    "LYN": {"letter": "K", "charge": 0, "tautomer": None},
    # Methionine
    "MET": {"letter": "M", "charge": 0, "tautomer": None},
    # Phenylalanine
    "PHE": {"letter": "F", "charge": 0, "tautomer": None},
    # Proline
    "PRO": {"letter": "P", "charge": 0, "tautomer": None},
    # Serine
    "SER": {"letter": "S", "charge": 0, "tautomer": None},
    # Threonine
    "THR": {"letter": "T", "charge": 0, "tautomer": None},
    # Tryptophan
    "TRP": {"letter": "W", "charge": 0, "tautomer": None},
    # Tyrosine
    "TYR": {"letter": "Y", "charge": 0, "tautomer": None},
    # Valine
    "VAL": {"letter": "V", "charge": 0, "tautomer": None},
}
