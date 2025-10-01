#!/usr/bin/env python3
"""
Helper script to get available methods from configuration.
Used by cli.sh for dynamic method validation.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from config.method_config import PYTHON_METHODS, R_METHODS

def main():
    """Print available methods."""
    if len(sys.argv) > 1 and sys.argv[1] == "--type":
        method_type = sys.argv[2] if len(sys.argv) > 2 else None
        if method_type == "python":
            print(" ".join(PYTHON_METHODS.keys()))
        elif method_type == "r":
            print(" ".join(R_METHODS.keys()))
        else:
            print("Error: Invalid method type", file=sys.stderr)
            sys.exit(1)
    else:
        # Print all methods
        all_methods = list(PYTHON_METHODS.keys()) + list(R_METHODS.keys())
        print(" ".join(all_methods))

if __name__ == "__main__":
    main()
