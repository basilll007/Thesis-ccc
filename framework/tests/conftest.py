"""Make the tests import the package from this source tree (not an installed copy)."""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
for m in [m for m in sys.modules if m == "ccc_audit" or m.startswith("ccc_audit.")]:
    del sys.modules[m]
