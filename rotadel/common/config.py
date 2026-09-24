import os
from pathlib import Path
import shutil

import pooch
import requests

ROOT_DIR = Path(__file__).resolve().parent.parent.parent
ANGLES_JSON_PATH = ROOT_DIR / "data" / "angles_combined.json"

# Config for xtb parallelisation
OMP_NUM_THREADS: int = 1
OMP_STACKSIZE: str = "1G"
OMP_MAX_ACTIVE_LEVELS: int = 1

# --- RotADeL database ---
# Two files are needed for RotADeL descriptors query:
# - the SQLite descriptors database (rotadel.db)
# - the NDRD file with backbone probabilities (ndrd_step10.csv)

# Get path to database directory with following priority:
# 1. Use ROTADEL_DATABASE_DIR environment variable if was set by user
if "ROTADEL_DATABASE_DIR" in os.environ:
    DATABASE_DIR = Path(os.environ["ROTADEL_DATABASE_DIR"])
# 2. Use database/ subfolder in the git repo if running from a clone
elif (ROOT_DIR / "pyproject.toml").exists():
    DATABASE_DIR = ROOT_DIR / "database"
# 3. Use user cache folder if `rotadel` installed via pip (pooch downloads database files there on first use, see below)
else:
    DATABASE_DIR = Path(pooch.os_cache("rotadel"))

SQL_NAME = "rotadel.db"
NDRD_NAME = "ndrd_step10.csv"
SQL_PATH = DATABASE_DIR / SQL_NAME
NDRD_PATH = DATABASE_DIR / NDRD_NAME

# RotADeL database files are stored on Zenodo, and downloaded automatically if not already present in DATABASE_DIR
_ZENODO_URL = "https://doi.org/10.5281/zenodo.22772020"


def _latest_zenodo_fetcher() -> pooch.Pooch:
    """Build a `pooch` fetcher for the latest published version of the database's Zenodo record."""
    latest_record_url = requests.head(_ZENODO_URL, allow_redirects=True, timeout=10).url
    record_id = latest_record_url.rstrip("/").rsplit("/", 1)[-1]
    metadata = requests.get(
        f"https://zenodo.org/api/records/{record_id}", timeout=10
    ).json()
    checksums = {f["key"]: f["checksum"] for f in metadata["files"]}
    return pooch.create(
        path=DATABASE_DIR,
        base_url=f"https://zenodo.org/records/{record_id}/files",
        registry=checksums,
    )


def fetch_database() -> None:
    """If database files are missing, download them from Zenodo."""
    # Skip if files already present
    if SQL_PATH.exists() and NDRD_PATH.exists():
        return
    print("\nFirst use - fetching the database from Zenodo...")
    try:
        fetcher = _latest_zenodo_fetcher()
        fetcher.fetch(SQL_NAME, progressbar=True)
        fetcher.fetch(NDRD_NAME, progressbar=True)
    except Exception as e:
        raise RuntimeError(
            f"Failed to fetch the database from Zenodo ({_ZENODO_URL}): {e}.\n"
            "See https://github.com/digital-chemistry-laboratory/rotadel#installation "
            "for how to download the database manually."
        ) from e
    print(f"Downloaded successfully to {DATABASE_DIR}")


def get_sql_path() -> Path:
    """Path to the SQL database file, fetching database from Zenodo on first use if not already present."""
    fetch_database()
    return SQL_PATH


def get_ndrd_path() -> Path:
    """Path to the NDRD csv file, fetching database from Zenodo on first use if not already present."""
    fetch_database()
    return NDRD_PATH


def clear_database_cache() -> None:
    """Delete rotadel's database cache folder."""
    cache_dir = Path(pooch.os_cache("rotadel"))
    if DATABASE_DIR != cache_dir:
        print(
            f"{DATABASE_DIR} is not rotadel's auto-managed cache folder, leaving it untouched."
        )
        return
    if cache_dir.exists():
        try:
            shutil.rmtree(cache_dir)
        except OSError:
            print(
                f"Could not fully remove {cache_dir}. Likely still used by another running process, "
                "e.g. a Jupyter kernel - try closing/restarting it, then run clear_database_cache() again."
            )
            return
        print(f"Removed {cache_dir}")


# ------
