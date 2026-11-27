"""Run or resume the final 40 nm shared-residual trajectory."""
from pathlib import Path
import sys


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
    from run_reset68_continuous import main

    main()
