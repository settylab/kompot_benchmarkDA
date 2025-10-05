#!/usr/bin/env python3
"""
Unified test runner for benchmarkDA.

Runs unit tests, integration tests, or all tests with proper environment setup.
"""

import sys
import subprocess
import argparse
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


def run_tests(test_type="all", verbose=False, pattern=None):
    """
    Run tests using pytest or fallback to direct execution.

    Parameters:
    -----------
    test_type : str
        Type of tests to run: "unit", "integration", or "all"
    verbose : bool
        Whether to show verbose output
    pattern : str
        Pattern to match test files (e.g., "test_critical*")

    Returns:
    --------
    int
        Exit code (0 = success, 1 = failure)
    """
    # Determine test directory
    test_dir = project_root / "tests"

    if test_type == "unit":
        test_path = test_dir / "unit"
    elif test_type == "integration":
        test_path = test_dir / "integration"
    else:  # all
        test_path = test_dir

    # Try to use pytest if available
    try:
        import pytest

        args = [str(test_path)]

        if verbose:
            args.append("-v")

        if pattern:
            args.extend(["-k", pattern])

        # Add color and better output
        args.extend(["--color=yes", "--tb=short"])

        print(f"\n{'=' * 70}")
        print(f"Running {test_type.upper()} tests using pytest")
        print(f"{'=' * 70}\n")

        return pytest.main(args)

    except ImportError:
        # Fallback to running tests directly
        print(f"\n{'=' * 70}")
        print(f"Running {test_type.upper()} tests (pytest not found, using direct execution)")
        print(f"{'=' * 70}\n")

        test_files = sorted(test_path.rglob("test_*.py"))

        if pattern:
            test_files = [f for f in test_files if pattern in f.name]

        if not test_files:
            print(f"No test files found in {test_path}")
            return 1

        failed = 0
        for test_file in test_files:
            print(f"\nRunning {test_file.name}...")
            print("-" * 70)

            result = subprocess.run(
                [sys.executable, str(test_file)],
                cwd=project_root,
            )

            if result.returncode != 0:
                failed += 1
                print(f"❌ {test_file.name} FAILED")
            else:
                print(f"✅ {test_file.name} PASSED")

        print(f"\n{'=' * 70}")
        if failed == 0:
            print(f"ALL TESTS PASSED ✅")
            print(f"{'=' * 70}\n")
            return 0
        else:
            print(f"{failed}/{len(test_files)} TEST(S) FAILED ❌")
            print(f"{'=' * 70}\n")
            return 1


def main():
    """Main entry point for test runner."""
    parser = argparse.ArgumentParser(
        description="Run benchmarkDA tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all tests
  python bin/run_tests.py

  # Run only unit tests
  python bin/run_tests.py --type unit

  # Run specific test pattern
  python bin/run_tests.py --pattern critical

  # Run with verbose output
  python bin/run_tests.py --verbose
        """,
    )

    parser.add_argument(
        "--type",
        choices=["unit", "integration", "all"],
        default="all",
        help="Type of tests to run (default: all)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show verbose output",
    )
    parser.add_argument(
        "-p",
        "--pattern",
        type=str,
        help="Pattern to match test files (e.g., 'critical')",
    )
    parser.add_argument(
        "--with-mamba",
        action="store_true",
        help="Run tests using mamba environment",
    )

    args = parser.parse_args()

    # If running with mamba, use mamba run
    if args.with_mamba:
        cmd = [
            "$MAMBA_EXE",
            "run",
            "-n",
            "benchmarkda",
            sys.executable,
            __file__,
            "--type",
            args.type,
        ]
        if args.verbose:
            cmd.append("--verbose")
        if args.pattern:
            cmd.extend(["--pattern", args.pattern])

        # Remove --with-mamba to avoid infinite recursion
        result = subprocess.run(" ".join(cmd), shell=True)
        sys.exit(result.returncode)

    # Run tests
    exit_code = run_tests(
        test_type=args.type, verbose=args.verbose, pattern=args.pattern
    )
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
