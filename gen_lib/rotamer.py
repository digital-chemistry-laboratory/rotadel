from typing import Any
from os import PathLike
from pathlib import Path
import json
import numpy as np

from gen_lib.optimisation import start_rotamer_xyz, get_opt_structures
from gen_lib.descriptors import calc_descriptors


class Rotamer:
    """Storing and manipulating data for a rotamer"""

    def __init__(
        self,
        key: str,
        dunbrack_data: dict[str, Any],
        charge: int,
        tautomer: str | None = None,
        run_path: str | PathLike | None = None,
    ) -> None:
        self._key = key
        self._dunbrack_data = dunbrack_data
        self._charge = charge
        self._tautomer = tautomer
        self._run_path = run_path if run_path is not None else Path.cwd()
        self._rotamer_elements = None
        self._rotamer_coordinates = None
        self._sidechainH_elements = None
        self._sidechainH_coordinates = None
        self._descriptors = {}

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
                same_chi_angles: bool = all(
                    self._dunbrack_data[field] == existing_rotamer[field]
                    for field in ["res", "chi1", "chi2", "chi3", "chi4"]
                )
                same_charge: bool = self._charge == existing_rotamer["charge"]
                same_tautomer: bool = self._tautomer == existing_rotamer["tautomer"]
                if same_chi_angles and same_charge and same_tautomer:
                    for name, entry in existing_rotamer.items():
                        if name not in self._dunbrack_data:
                            if name == "rotamer":
                                self._rotamer_elements = np.array(entry["elements"])
                                self._rotamer_coordinates = np.array(
                                    entry["coordinates"]
                                )
                            elif name == "sidechain-H":
                                self._sidechainH_elements = np.array(entry["elements"])
                                self._sidechainH_coordinates = np.array(
                                    entry["coordinates"]
                                )
                            else:
                                self._descriptors[name] = entry
                    return True
        return False

    def opt_geometries(self) -> None:
        """Optimise the rotamer and sidechain-H geometries."""
        starting_rotamer_xyz = self._run_path / "rotamer_start.xyz"
        start_rotamer_xyz(
            dunbrack_data=self._dunbrack_data,
            xyz_output=starting_rotamer_xyz,
            charge=self._charge,
            tautomer=self._tautomer,
        )
        el_whole, coord_whole, el_sidechain, coord_sidechain = get_opt_structures(
            dunbrack_data=self._dunbrack_data,
            charge=self._charge,
            run_folder=self._run_path,
            rotamer_xyz=starting_rotamer_xyz,
        )
        self._rotamer_elements = el_whole
        self._rotamer_coordinates = coord_whole
        self._sidechainH_elements = el_sidechain
        self._sidechainH_coordinates = coord_sidechain

    def calc_descriptors(self) -> None:
        """Calculate the descriptors on the sidechain-H."""
        if self._sidechainH_elements is None or self._sidechainH_coordinates is None:
            raise ValueError(
                "The geometry must first be optimsed before calculating descriptors."
            )
        is_pro = self._dunbrack_data["res"] == "PRO"
        descriptors = calc_descriptors(
            self._sidechainH_elements,
            self._sidechainH_coordinates,
            charge=self._charge,
            is_pro=is_pro,
        )
        self._descriptors = descriptors

    def to_dict(self) -> dict[str, Any]:
        """Return a json-compatible dictionary with the rotamer data."""
        dictionary = {
            self._key: {
                **self._dunbrack_data,
                "charge": self._charge,
                "tautomer": self._tautomer,
                "rotamer": {
                    "elements": self._rotamer_elements.tolist(),
                    "coordinates": self._rotamer_coordinates.tolist(),
                },
                "sidechain-H": {
                    "elements": self._sidechainH_elements.tolist(),
                    "coordinates": self._sidechainH_coordinates.tolist(),
                },
            },
        }
        if self._descriptors is not None:
            for descriptor_name, descriptor_value in self._descriptors.items():
                dictionary[self._key][descriptor_name] = descriptor_value
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
