"""Locate the checkout containing user presets, caches, and generated maps."""
import os
from pathlib import Path


def project_root():
    explicit = os.environ.get("NZ_TERRAIN_ROOT")
    if explicit:
        root = Path(explicit).expanduser().resolve()
        if not (root / "locations.json").is_file():
            raise FileNotFoundError(f"NZ_TERRAIN_ROOT has no locations.json: {root}")
        return root
    checkout = Path(__file__).resolve().parents[2]
    for root in [checkout, Path.cwd(), *Path.cwd().parents]:
        if (root / "locations.json").is_file() and (root / "water-configs").is_dir():
            return root
    raise FileNotFoundError("Run from an nz-terrain checkout or set NZ_TERRAIN_ROOT to its path")


def preview_font_path():
    """Use the dependency's bundled font instead of an OS-specific font path."""
    import matplotlib
    return str(Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans.ttf")
