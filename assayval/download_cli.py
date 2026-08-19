"""Console entry point for the assembly-download workflow.

Thin wrapper that forwards all arguments to the bundled
scripts/download/download_assemblies.sh. That script self-locates its helper
parse_metadata.py and the optional credentials file (via BASH_SOURCE), so a
single command runs the whole download-and-parse-metadata workflow from any
working directory. Exposed as the `download-assemblies` command.
"""
import subprocess
import sys
from pathlib import Path

# <repo>/assayval/download_cli.py -> <repo>/scripts/download/download_assemblies.sh
SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "download" / "download_assemblies.sh"


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if not SCRIPT.exists():
        print(f"error: download script not found at {SCRIPT}", file=sys.stderr)
        return 1
    return subprocess.run(["bash", str(SCRIPT), *args]).returncode


if __name__ == "__main__":
    sys.exit(main())
