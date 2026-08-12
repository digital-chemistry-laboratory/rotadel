"""Extract all rotamers of one amino acid / one backbone conformation from the database
as backbone-aligned xyz files (one for the backbone, one sidechain per rotamer).

Run: python extract_ensemble.py, then python chimerax_script.py
"""

from pathlib import Path
import sqlite3
import numpy as np

from rotadel.common.config import SQL_PATH

RES = "MET"
OUT_DIR = Path(__file__).parent / f"data_{RES}"
# MET: 3 chi angles (enough variation), only 4 sidechain heavy atoms (readable overlay),
# S atom gives a visual anchor along the chain -> 27 rotamers for a given backbone
CHARGE = 0
TAUTOMER = None
# fully extended backbone (beta region) -> classical zigzag
PHI, PSI = -180.0, -180.0
MIN_PROB = 0.0  # keep all rotamers of the bin
WRITE_H = False  # heavy atoms only, much cleaner overlay


def write_xyz(path, elements, coords, keep, comment=None):
    """Write the atoms of indices `keep` to an xyz file (H skipped if WRITE_H is False)."""
    keep = [i for i in keep if WRITE_H or elements[i] != "H"]
    lines = [f"{len(keep)}", comment or path.stem]
    lines += [
        f"{elements[i]:2s} {' '.join(f'{v:12.6f}' for v in coords[i])}" for i in keep
    ]
    path.write_text("\n".join(lines) + "\n")


def extract_ensemble():
    """Write the backbone and one backbone-aligned sidechain xyz file per rotamer."""
    con = sqlite3.connect(SQL_PATH)
    ids = con.execute(
        "SELECT rotamer_id, prob FROM rotamers_data WHERE res=? AND phi=? AND psi=? "
        "AND charge=? AND tautomer IS ? AND prob>=? ORDER BY prob DESC",
        (RES, PHI, PSI, CHARGE, TAUTOMER, MIN_PROB),
    ).fetchall()
    if not ids:
        raise SystemExit(f"no rotamer found for {RES} phi={PHI} psi={PSI}")

    def geometry(rotamer_id):
        """Return (elements, coordinates) of a whole rotamer, in database atom order."""
        rows = con.execute(
            "SELECT element, x, y, z FROM rotamer_xyz WHERE rotamer_id=? ORDER BY atom_idx",
            (rotamer_id,),
        ).fetchall()
        return [r[0] for r in rows], np.array([r[1:] for r in rows], dtype=float)

    # atom order is always N,H,H,CA,HA,C,O | sidechain | OXT,HXT
    ref_elements, ref_coords = geometry(ids[0][0])
    n = len(ref_elements)
    bb_idx = [0, 1, 2, 3, 4, 5, 6, n - 2, n - 1]
    align_idx = [0, 3, 4, 5]  # N, CA, HA, C
    sc_idx = [3] + list(range(7, n - 2))  # CA kept so that the CA-CB bond is drawn

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_xyz(OUT_DIR / "backbone.xyz", ref_elements, ref_coords, bb_idx)
    rmsds = []
    # rank in the file name so that sorting them = decreasing probability,
    # probability in the xyz comment line so that the ChimeraX script can weight them
    for rank, (rotamer_id, prob) in enumerate(ids):
        elements, coords = geometry(rotamer_id)
        # Kabsch superposition of the backbone onto the most probable rotamer
        p, q = coords[align_idx], ref_coords[align_idx]
        pc, qc = p.mean(0), q.mean(0)
        u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
        rot = vt.T @ np.diag([1, 1, np.sign(np.linalg.det(vt.T @ u.T))]) @ u.T
        coords = (coords - pc) @ rot.T + qc
        rmsds.append(np.sqrt(((coords[align_idx] - q) ** 2).sum(1).mean()))
        write_xyz(
            OUT_DIR / f"sc_{rank:02d}_{rotamer_id}.xyz",
            elements,
            coords,
            sc_idx,
            comment=f"{rotamer_id} prob={prob:.6f}",
        )

    con.close()
    print(
        f"{len(ids)} rotamers -> {OUT_DIR}, backbone alignment RMSD max {max(rmsds):.3f} A"
    )


if __name__ == "__main__":
    extract_ensemble()
