from typing import Any
from os import PathLike
from pathlib import Path
import json
import sqlite3
import numpy as np

from aa_descriptors_library.optimisation import start_rotamer_xyz, get_opt_structures
from aa_descriptors_library.descriptors import calc_descriptors


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

    def load_existing_sidechain_json(self, json_file: str | PathLike) -> bool:
        """Check if the sidechain already exists in the JSON file, if yes load its coordinates and descriptors.
        Args:
            JSON file with already calculated rotamer data
        Returns:
            If same sidechain is already in JSON:
                True and copy existing sidechain and descriptors in Rotamer instance attributes.
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
                        if name == "sidechain-H":
                            self._sidechainH_elements = np.array(entry["elements"])
                            self._sidechainH_coordinates = np.array(
                                entry["coordinates"]
                            )
                        elif name == "descriptors":
                            self._descriptors.update(entry)
                    return True
        return False

    def load_existing_sidechain_sql(self, conn: sqlite3.Connection) -> bool:
        """If a matching sidechain already exists in the SQLite DB, if yes load its coordinates and descriptors.
        Args:
            conn: SQLite connection to the database
        Returns:
            If same sidechain is already in DB:
                True and copy existing sidechain and descriptors in Rotamer instance attributes.
            If not: False and does not copy data.
        """
        cur = conn.cursor()
        cur.execute(
            """
            SELECT rotamer_id, descriptors
            FROM rotamers_data
            WHERE res = ? AND chi2 = ? AND chi3 IS ? AND chi4 IS ?
            AND charge = ? AND tautomer IS ?
            LIMIT 1
        """,
            (
                self._dunbrack_data["res"],
                self._dunbrack_data.get("chi2"),
                self._dunbrack_data.get("chi3"),
                self._dunbrack_data.get("chi4"),
                self._charge,
                self._tautomer,
            ),
        )
        match = cur.fetchone()
        if not match:
            return False

        matching_rotamer_key, descriptors = match

        cur.execute(
            """
            SELECT element, x, y, z
            FROM sidechainH_xyz
            WHERE rotamer_id = ?
            ORDER BY atom_idx
        """,
            (matching_rotamer_key,),
        )
        rows = cur.fetchall()
        self._sidechainH_elements = np.array([r[0] for r in rows])
        self._sidechainH_coordinates = np.array([r[1:] for r in rows], dtype=float)

        self._descriptors.update(json.loads(descriptors))

        return True

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
        descriptors = calc_descriptors(
            self._sidechainH_elements,
            self._sidechainH_coordinates,
            charge=self._charge,
            res=self._dunbrack_data["res"],
        )
        self._descriptors = descriptors

    def to_dict(self) -> dict[str, Any]:
        """Return a json-compatible dictionary with the rotamer data."""
        dictionary = {
            self._key: {
                **self._dunbrack_data,
                "charge": self._charge,
                "tautomer": self._tautomer,
            },
        }
        if self._rotamer_coordinates is not None:
            dictionary[self._key]["rotamer"] = {
                "elements": self._rotamer_elements.tolist(),
                "coordinates": self._rotamer_coordinates.tolist(),
            }
        dictionary[self._key]["sidechain-H"] = {
            "elements": self._sidechainH_elements.tolist(),
            "coordinates": self._sidechainH_coordinates.tolist(),
        }
        if self._descriptors:
            dictionary[self._key]["descriptors"] = {}
            for descriptor_name, descriptor_value in self._descriptors.items():
                dictionary[self._key]["descriptors"][descriptor_name] = descriptor_value
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

    def to_sql(self, conn: sqlite3.Connection) -> None:
        """Save the rotamer data in a SQLite database."""
        cur = conn.cursor()
        cur.execute(
            """
            INSERT OR REPLACE INTO rotamers_data (rotamer_id, res, letter, phi, psi,
            chi1, chi2, chi3, chi4, prob, charge, tautomer, descriptors)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
            (
                self._key,
                self._dunbrack_data["res"],
                self._dunbrack_data["letter"],
                self._dunbrack_data["phi"],
                self._dunbrack_data["psi"],
                self._dunbrack_data["chi1"],
                self._dunbrack_data.get("chi2"),
                self._dunbrack_data.get("chi3"),
                self._dunbrack_data.get("chi4"),
                self._dunbrack_data["prob"],
                self._charge,
                self._tautomer,
                json.dumps(self._descriptors),
            ),
        )

        if self._rotamer_coordinates is not None:
            for i, (el, coord) in enumerate(
                zip(self._rotamer_elements, self._rotamer_coordinates)
            ):
                cur.execute(
                    """
                    INSERT INTO rotamer_xyz (rotamer_id, atom_idx, element, x, y, z)
                    VALUES (?, ?, ?, ?, ?, ?)
                """,
                    (self._key, i, el, *coord),
                )

        for i, (el, coord) in enumerate(
            zip(self._sidechainH_elements, self._sidechainH_coordinates)
        ):
            cur.execute(
                """
                INSERT INTO sidechainH_xyz (rotamer_id, atom_idx, element, x, y, z)
                VALUES (?, ?, ?, ?, ?, ?)
            """,
                (self._key, i, el, *coord),
            )

        conn.commit()
