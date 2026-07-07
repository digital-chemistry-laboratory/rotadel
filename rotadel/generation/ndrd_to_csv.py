"""
Save the probability data from NDRD library into a CSV file.
"""

import pandas as pd

from rotadel.common.config import NDRD_PATH, ROOT_DIR

lib_file = ROOT_DIR / "data" / "Dunbrack_libraries" / "ndrd" / "NDRD_TCBIG.txt"

column_names = [
    "res",
    "neighbour_pos",
    "neighbour_res",
    "phi",
    "psi",
    "prob",
    "log_prob",
    "cum_sum",
]
lib_df = pd.read_csv(
    lib_file,
    sep=r"\s+",
    comment="#",
    names=column_names,
    low_memory=False,
)
lib_df.drop(["cum_sum"], axis="columns", inplace=True)

# Keep only the phi/psi incremented by 10° (increment of 10° in rotamer library while 5° in NDRD)
lib_df_step_10 = lib_df.loc[
    (lib_df["phi"] % 10 == 0) & (lib_df["psi"] % 10 == 0)
].copy()

lib_df_step_10.to_csv(NDRD_PATH, index=False)
