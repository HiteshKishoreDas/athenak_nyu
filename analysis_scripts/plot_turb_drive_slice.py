#!/usr/bin/env python3

import os
import sys
import gc
import numpy as np
import cmasher as cm
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib"))

import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import bin_convert  # noqa: E402

# # Edit these directly or override them before calling main().
# SNAPSHOT: Path | None = None
# RUN_DIR: Path = REPO_ROOT / "turb_drive"
# VAR: str = "dens"
# AXIS: str = "z"
# INDEX: str | int = "mid"
# OUTPUT: Path | None = None
# LINEAR: bool = False
# CMAP: str = "magma"


def load_snapshot(path: Path) -> dict[str, object]:
    return bin_convert.read_binary_as_athdf(str(path))


snapshot = "bin/Turb.slice_x2.00040.bin"
data = load_snapshot(snapshot)

plt.figure()
# plt.imshow(np.log10(data["dens"])[:, 0, :], vmin=-1, vmax=1)
plt.imshow(data["velz"][:, 0, :], cmap=cm.redshift_r, vmin=-1, vmax=1)
plt.colorbar()
plt.show()

del data
gc.collect()
