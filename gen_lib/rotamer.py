from typing import Any
from os import PathLike
from pathlib import Path
import json
import numpy as np

from gen_lib.optimisation import init_rotamer_xyz, get_opt_structures
from gen_lib.descriptors import calc_descriptors


class Rotamer:
    """Storing and manipulating data for a rotamer"""

    def __init__(
        self,
        key: str,
        dunbrack_data: dict[str, Any] | None = None,
        charge: int | None = None,
        run_path: str | PathLike | None = None,
    ) -> None:
        self.key = key
        self.dunbrack_data = dunbrack_data
        self.charge = charge
        self.run_path = run_path if run_path is not None else Path.cwd()
        self.rotamer_elements = None
        self.rotamer_coordinates = None
        self.sidechainH_elements = None
        self.sidechainH_coordinates = None
        self.descriptors = {}

    def load_existing_sidechain(self, json_file: str | PathLike) -> bool:
        """Check if the sidechain already exists in the JSON file and load calculated data if applicable.
        Args:
            JSON file with already calculated rotamer data
        Returns:
            If same sidechain is already in JSON: True and copy existing data in Rotamer instance attributes.
            If not: False and does not copy data.
        """
        if Path(json_file).exists():
            with open(json_file, "r") as f:
                content = json.load(f)
            for _, existing_rotamer in content.items():
                if all(
                    self.dunbrack_data[field] == existing_rotamer[field]
                    for field in ["res", "chi1", "chi2", "chi3", "chi4"]
                ):
                    for name, entry in existing_rotamer.items():
                        if name not in self.dunbrack_data:
                            if name == "rotamer":
                                self.rotamer_elements = np.array(entry["elements"])
                                self.rotamer_coordinates = np.array(
                                    entry["coordinates"]
                                )
                            elif name == "sidechain-H":
                                self.sidechainH_elements = np.array(entry["elements"])
                                self.sidechainH_coordinates = np.array(
                                    entry["coordinates"]
                                )
                            else:
                                self.descriptors[name] = entry
                    return True
        return False

    def opt_geometries(self) -> None:
        """Optimise the rotamer and sidechain-H geometries."""
        starting_rotamer_xyz = self.run_path / "rotamer_start.xyz"
        init_rotamer_xyz(
            rotamer_id=self.key,
            dunbrack_data=self.dunbrack_data,
            xyz_output=starting_rotamer_xyz,
        )
        el_whole, coord_whole, el_sidechain, coord_sidechain = get_opt_structures(
            dunbrack_data=self.dunbrack_data,
            charge=self.charge,
            run_folder=self.run_path,
            rotamer_xyz=starting_rotamer_xyz,
        )
        self.rotamer_elements = el_whole
        self.rotamer_coordinates = coord_whole
        self.sidechainH_elements = el_sidechain
        self.sidechainH_coordinates = coord_sidechain

    def calc_descriptors(self) -> None:
        """Calculate the descriptors on the sidechain-H."""
        if self.sidechainH_elements is None or self.sidechainH_coordinates is None:
            raise ValueError(
                "The geometry must first be optimsed before calculating descriptors."
            )
        descriptors = calc_descriptors(
            self.sidechainH_elements, self.sidechainH_coordinates, charge=self.charge
        )
        self.descriptors = descriptors

    def to_dict(self) -> dict[str, Any]:
        """Return a json-compatible dictionary with the rotamer data."""
        dictionary = {
            self.key: {
                **self.dunbrack_data,
                "rotamer": {
                    "elements": self.rotamer_elements.tolist(),
                    "coordinates": self.rotamer_coordinates.tolist(),
                },
                "sidechain-H": {
                    "elements": self.sidechainH_elements.tolist(),
                    "coordinates": self.sidechainH_coordinates.tolist(),
                },
            },
        }
        if self.descriptors is not None:
            for descriptor_name, descriptor_value in self.descriptors.items():
                dictionary[self.key][descriptor_name] = descriptor_value
        return dictionary

    def to_json(self, json_file: str | PathLike) -> None:
        """Save the rotamer data in a JSON file."""
        if Path(json_file).exists():
            with open(json_file, "r") as f:
                data = json.load(f)
        else:
            data = {}
        data.update(self.to_dict())
        with open(json_file, "w") as f:
            json.dump(data, f, indent=4)
