from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent.parent
SQL_PATH = PACKAGE_DIR / "database" / "rotadel.db"
NDRD_PATH = PACKAGE_DIR / "database" / "ndrd_step10.csv"
ROOT_DIR = Path(__file__).resolve().parent.parent.parent
ANGLES_JSON_PATH = ROOT_DIR / "data" / "angles_combined.json"

# Config for xtb parallelisation
OMP_NUM_THREADS: int = 1
OMP_STACKSIZE: str = "1G"
OMP_MAX_ACTIVE_LEVELS: int = 1
