class Rotamer:
    """Storing and manipulating data for a rotamer"""

    def __init__(self, key: str):
        self.key = key
        self.rotamer_elements = None
        self.rotamer_coordinates = None
        self.sidechainH_elements = None
        self.sidechainH_coordinates = None

    def set_opt_geometries(
        self,
        rotamer_elements: list,
        rotamer_coordinates: list,
        sidechainH_elements: list,
        sidechainH_coordinates: list,
    ) -> None:
        """Set the elements and coordinates of the optimised rotamer and sidechain-H geometries."""
        self.rotamer_elements = rotamer_elements
        self.rotamer_coordinates = rotamer_coordinates
        self.sidechainH_elements = sidechainH_elements
        self.sidechainH_coordinates = sidechainH_coordinates

    def to_dict(self):
        """Return a dictionary with the rotamer data."""
        return {
            self.key: {
                "rotamer": {
                    "elements": self.rotamer_elements,
                    "coordinates": self.rotamer_coordinates,
                },
                "sidechain-H": {
                    "elements": self.sidechainH_elements,
                    "coordinates": self.sidechainH_coordinates,
                },
            }
        }
