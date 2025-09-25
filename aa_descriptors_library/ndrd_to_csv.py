"""
Save the probability data from NDRD library into a CSV file.
"""

import pandas as pd

lib_file = "/cluster/project/jorner/lajacot/projects/aa-descriptors-library/Dunbrack_library/ndrd/NDRD_TCBIG.txt"
csv_output = "/cluster/project/jorner/lajacot/projects/aa-descriptors-library/gen_lib/data/ndrd.csv"

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
lib_df.to_csv(csv_output, index=False)
