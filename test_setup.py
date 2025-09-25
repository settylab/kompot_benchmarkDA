#!/usr/bin/env python3
"""
Comprehensive test script for BenchmarkDA setup.
Tests environment, dependencies, and basic functionality.
"""

import sys
import subprocess
import importlib
import os
from pathlib import Path


def print_colored(text, color):
    colors = {
        'red': '\033[0;31m',
        'green': '\033[0;32m',
        'yellow': '\033[1;33m',
        'blue': '\033[0;34m',
        'bold': '\033[1m',
        'reset': '\033[0m'
    }
    print(f"{colors.get(color, '')}{text}{colors['reset']}")

def test_python_packages():
    """Test Python package imports."""
    print_colored("Testing Python packages...", 'blue')

    required_packages = {
        'numpy': 'numpy',
        'pandas': 'pandas',
        'scanpy': 'scanpy',
        'anndata': 'anndata',
        'sklearn': 'scikit-learn',
        'matplotlib': 'matplotlib',
        'meld': 'meld',
        'cna': 'cna',
        'palantir': 'palantir',
        'mellon': 'mellon',
        'kompot': 'kompot'
    }

    success_count = 0
    for package, display_name in required_packages.items():
        try:
            importlib.import_module(package)
            print_colored(f"✓ {display_name}", 'green')
            success_count += 1
        except ImportError as e:
            print_colored(f"✗ {display_name}: {e}", 'red')

    print_colored(f"\nPython packages: {success_count}/{len(required_packages)} available",
                  'green' if success_count == len(required_packages) else 'yellow')
    return success_count == len(required_packages)

def test_r_packages():
    """Test R package availability."""
    print_colored("\nTesting R packages...", 'blue')

    r_packages = ['argparse', 'tidyverse', 'SingleCellExperiment', 'scran',
                  'Seurat', 'igraph']

    cmd = ['R', '--vanilla', '--slave', '-e', f'''
packages <- c({", ".join(f'"{pkg}"' for pkg in r_packages)})
missing <- c()
for (pkg in packages) {{
    if (!requireNamespace(pkg, quietly = TRUE)) {{
        missing <- c(missing, pkg)
        cat(paste("✗", pkg, "\\n"))
    }} else {{
        cat(paste("✓", pkg, "\\n"))
    }}
}}
if (length(missing) > 0) {{
    cat(paste("\\nMissing R packages:", paste(missing, collapse = ", "), "\\n"))
}} else {{
    cat("\\nAll R packages available!\\n")
}}
    ''']

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print_colored(f"R warnings: {result.stderr}", 'yellow')
        return "Missing R packages:" not in result.stdout
    except Exception as e:
        print_colored(f"✗ R not available: {e}", 'red')
        return False

def test_environment_detection():
    """Test environment detection utilities."""
    print_colored("\nTesting environment detection...", 'blue')

    try:
        result = subprocess.run(['bash', 'bin/environment_utils.sh'],
                              capture_output=True, text=True,
                              cwd=Path(__file__).parent)

        # Test individual functions
        test_cmd = '''
        source bin/environment_utils.sh
        echo "Package manager: $(detect_package_manager)"
        echo "Environment: $(detect_benchmarkda_environment)"
        '''

        result = subprocess.run(['bash', '-c', test_cmd],
                              capture_output=True, text=True,
                              cwd=Path(__file__).parent)

        print(result.stdout)
        return "Package manager:" in result.stdout

    except Exception as e:
        print_colored(f"✗ Environment detection failed: {e}", 'red')
        return False

def test_file_structure():
    """Test that all required files exist."""
    print_colored("\nTesting file structure...", 'blue')

    required_files = [
        'main.sh',
        'setup_environment.sh',
        'bin/environment_utils.sh',
        'bin/run_benchmark.py',
        'run_da_method.py',
        'config/method_config.py',
        'config/dataset_config.py',
        'python_method/shared_embedding_utils.py',
        'python_method/data_loader.py',
        'scripts/run_DA.r',
        'environment_minimal.yml',
        'environment_complete.yml'
    ]

    project_root = Path(__file__).parent
    missing_files = []

    for file_path in required_files:
        full_path = project_root / file_path
        if full_path.exists():
            print_colored(f"✓ {file_path}", 'green')
        else:
            print_colored(f"✗ {file_path}", 'red')
            missing_files.append(file_path)

    if missing_files:
        print_colored(f"\nMissing files: {', '.join(missing_files)}", 'red')
        return False
    else:
        print_colored("\nAll required files present!", 'green')
        return True

def test_method_configs():
    """Test method configuration validity."""
    print_colored("\nTesting method configurations...", 'blue')

    try:
        sys.path.append(str(Path(__file__).parent))
        from config.method_config import PYTHON_METHODS, R_METHODS

        print_colored(f"✓ Python methods: {', '.join(PYTHON_METHODS.keys())}", 'green')
        print_colored(f"✓ R methods: {', '.join(R_METHODS.keys())}", 'green')

        # Check that method scripts exist
        project_root = Path(__file__).parent
        for method, config in PYTHON_METHODS.items():
            script_path = project_root / 'python_method' / config['script']
            if script_path.exists():
                print_colored(f"✓ {config['script']}", 'green')
            else:
                print_colored(f"✗ {config['script']}", 'red')
                return False

        return True

    except Exception as e:
        print_colored(f"✗ Method config error: {e}", 'red')
        return False

def run_basic_functionality_test():
    """Test basic functionality with a simple method."""
    print_colored("\nTesting basic functionality...", 'blue')

    # Test environment detection
    try:
        # Use MAMBA_EXE to run with correct environment
        mamba_exe = os.environ.get('MAMBA_EXE', 'mamba')
        result = subprocess.run([mamba_exe, 'run', '-n', 'benchmarkda', 'python', 'run_da_method.py', '--list'],
                              capture_output=True, text=True,
                              cwd=Path(__file__).parent)

        if "AVAILABLE METHODS" in result.stdout:
            print_colored("✓ Method listing works", 'green')
            return True
        else:
            print_colored(f"✗ Method listing failed: {result.stderr}", 'red')
            return False

    except Exception as e:
        print_colored(f"✗ Basic functionality test failed: {e}", 'red')
        return False

def main():
    """Run all tests."""
    print_colored("=" * 60, 'bold')
    print_colored("BenchmarkDA Setup Test Suite", 'bold')
    print_colored("=" * 60, 'bold')

    tests = [
        ("File Structure", test_file_structure),
        ("Method Configs", test_method_configs),
        ("Environment Detection", test_environment_detection),
        ("Python Packages", test_python_packages),
        ("R Packages", test_r_packages),
        ("Basic Functionality", run_basic_functionality_test)
    ]

    results = {}
    for test_name, test_func in tests:
        print_colored(f"\n{'=' * 40}", 'blue')
        print_colored(f"Running: {test_name}", 'bold')
        print_colored(f"{'=' * 40}", 'blue')

        try:
            results[test_name] = test_func()
        except Exception as e:
            print_colored(f"✗ {test_name} crashed: {e}", 'red')
            results[test_name] = False

    # Summary
    print_colored("\n" + "=" * 60, 'bold')
    print_colored("TEST SUMMARY", 'bold')
    print_colored("=" * 60, 'bold')

    passed = sum(results.values())
    total = len(results)

    for test_name, result in results.items():
        status = "PASS" if result else "FAIL"
        color = 'green' if result else 'red'
        print_colored(f"{test_name:<25} {status}", color)

    print_colored(f"\nOverall: {passed}/{total} tests passed",
                  'green' if passed == total else 'yellow')

    if passed == total:
        print_colored("🎉 All tests passed! BenchmarkDA is ready to use.", 'green')
        print_colored("\nNext steps:", 'blue')
        print_colored("1. Run: bash main.sh", 'blue')
        print_colored("2. Or run individual methods: python run_da_method.py --list", 'blue')
    else:
        print_colored("⚠️  Some tests failed. Please fix issues before proceeding.", 'yellow')
        if not results.get("Python Packages", True):
            print_colored("Run: bash setup_environment.sh", 'yellow')

    return passed == total

if __name__ == '__main__':
    # Check if we're in a reasonable Python environment
    if sys.version_info < (3, 6):
        print("ERROR: This script requires Python 3.6+")
        print("To run properly, use: mamba run -n benchmarkda python test_setup.py")
        print("Or: $MAMBA_EXE run -n benchmarkda python test_setup.py")
        sys.exit(1)

    success = main()
    sys.exit(0 if success else 1)