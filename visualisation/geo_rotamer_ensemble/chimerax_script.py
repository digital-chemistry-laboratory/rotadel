"""Write the ChimeraX script overlaying the sidechains of a rotamer ensemble.

Run: python chimerax_script.py [ensemble_folder]
Then open the generated view_ensemble.cxc in ChimeraX
and save image with: `save path/image.png supersample 4`
"""

from pathlib import Path
import sys

RES = "MET"
ENSEMBLE_DIR = Path(__file__).parent / f"{RES}_ext"

# Transparency (%) is spread over rotamer ranks
MIN_TRANSPARENCY = 0
# transparency of the second most probable rotamer => isolates the first rotamer visually
SECOND_TRANSPARENCY = 65
# transparency of the least probable rotamer
MAX_TRANSPARENCY = 95
# fade of the rotamers after the first one, < 1 steep at the start, 1 is linear
GAMMA = 0.7

# carbons of the backbone
BB_COLOR = "gray"
# carbons of the sidechains
SC_COLOR = "tan"


def write_chimerax_script(ensemble_dir):
    """Write view_ensemble.cxc in `ensemble_dir`, opening backbone as #1 and sidechains after."""
    sc_files = sorted(ensemble_dir.glob("sc_*.xyz"))  # sorted by rank = decreasing prob
    last = len(sc_files) + 1
    span, n = MAX_TRANSPARENCY - SECOND_TRANSPARENCY, max(len(sc_files) - 2, 1)
    fade = [MIN_TRANSPARENCY] + [
        SECOND_TRANSPARENCY + span * ((i - 1) / n) ** GAMMA
        for i in range(1, len(sc_files))
    ]
    cxc = ["close session", "open backbone.xyz"]
    cxc += [f"open {f.name}" for f in sc_files]
    cxc += [
        "style ~H ball",
        f"color #1 {BB_COLOR}",
        f"color #2-{last} {SC_COLOR}",
        # byhetero restores the N/O/S element colours
        f"color #1-{last} byhetero",
        # the CA is duplicated in every sidechain file but belongs to the backbone
        # (C1 = first carbon of the file, ChimeraX names xyz atoms element by element)
        f"color #2-{last}@C1 {BB_COLOR}",
    ]
    # "target ab" is mandatory: the transparency command defaults to surfaces only
    cxc += [f"transparency #{i + 2} {t:.0f} target ab" for i, t in enumerate(fade)]
    cxc += ["set bgColor white", "view"]
    script = ensemble_dir / "view_ensemble.cxc"
    script.write_text("\n".join(cxc) + "\n")
    print(f"{script} ({last} models)")


if __name__ == "__main__":
    write_chimerax_script(Path(sys.argv[1]) if len(sys.argv) > 1 else ENSEMBLE_DIR)
