from typing import Any
from morfeus.typing import Array1DStr, Array2DFloat


class Rotamer:
    """Storing and manipulating data for a rotamer"""

    def __init__(self, key: str) -> None:
        self.key = key
        self.rotamer_elements = None
        self.rotamer_coordinates = None
        self.sidechainH_elements = None
        self.sidechainH_coordinates = None
        self.descriptors = None

    def set_opt_geometries(
        self,
        rotamer_elements: Array1DStr,
        rotamer_coordinates: Array2DFloat,
        sidechainH_elements: Array1DStr,
        sidechainH_coordinates: Array2DFloat,
    ) -> None:
        """Set the elements and coordinates of the optimised rotamer and sidechain-H geometries."""
        self.rotamer_elements = rotamer_elements
        self.rotamer_coordinates = rotamer_coordinates
        self.sidechainH_elements = sidechainH_elements
        self.sidechainH_coordinates = sidechainH_coordinates

    def to_dict(self) -> dict[str, Any]:
        """Return a json-compatible dictionary with the rotamer data."""
        dictionary = {
            self.key: {
                "rotamer": {
                    "elements": self.rotamer_elements.tolist(),
                    "coordinates": self.rotamer_coordinates.tolist(),
                },
                "sidechain-H": {
                    "elements": self.sidechainH_elements.tolist(),
                    "coordinates": self.sidechainH_coordinates.tolist(),
                },
            }
        }
        if self.descriptors is not None:
            for descriptor_name, descriptor_value in self.descriptors.items():
                dictionary[self.key][descriptor_name] = descriptor_value
        return dictionary
