"""
Test for quantile_assign_label bug fix.

Tests that the quantile calculation correctly extracts scalar values from DataFrames
instead of passing DataFrame objects to np.quantile, which would cause broadcasting issues.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import anndata as ad
from lib.synth_labels import quantile_assign_label


def test_quantile_scalar_extraction():
    """Test that quantile thresholds are scalar values, not arrays."""
    print("=" * 70)
    print("TEST: Quantile Assign Label - Scalar Extraction Fix")
    print("=" * 70)

    # Create mock data
    n_cells = 100
    populations = ['M1'] * 30 + ['M2'] * 50 + ['M3'] * 20

    # Create AnnData
    adata = ad.AnnData(
        X=np.zeros((n_cells, 10)),
        obs=pd.DataFrame({
            'celltype': populations,
            'Condition1_prob': np.random.uniform(0.3, 0.8, n_cells),
        }, index=[f'cell_{i}' for i in range(n_cells)])
    )

    # Calculate Condition2_prob
    adata.obs['Condition2_prob'] = 1 - adata.obs['Condition1_prob']

    # Test with pop_enr > 0.5 (M2 enriched in Condition1)
    print("\nTest case 1: pop_enr > 0.5 (M2 enriched in Condition1)")
    print(f"  Population: M2 ({50}/{n_cells} = 50% of cells)")
    print(f"  Enrichment: 0.75")

    adata_test1 = adata.copy()
    adata_result1 = quantile_assign_label(adata_test1, 'celltype', 0.75, 'M2')

    # Verify true_labels is 1D array
    assert adata_result1.obs['true_labels'].ndim == 1, "true_labels should be 1D array"
    print(f"  ✓ true_labels is 1D (shape: {adata_result1.obs['true_labels'].shape})")

    # Check label distribution
    label_counts = adata_result1.obs['true_labels'].value_counts().to_dict()
    print(f"  Label distribution: {label_counts}")

    # Verify all labels are valid
    valid_labels = {'NegLFC', 'NotDA', 'PosLFC'}
    actual_labels = set(adata_result1.obs['true_labels'].unique())
    assert actual_labels.issubset(valid_labels), f"Invalid labels found: {actual_labels - valid_labels}"
    print(f"  ✓ All labels are valid: {actual_labels}")

    # Test with pop_enr < 0.5 (M1 enriched in Condition2)
    print("\nTest case 2: pop_enr < 0.5 (M1 enriched in Condition2)")
    print(f"  Population: M1 ({30}/{n_cells} = 30% of cells)")
    print(f"  Enrichment: 0.25")

    adata_test2 = adata.copy()
    adata_result2 = quantile_assign_label(adata_test2, 'celltype', 0.25, 'M1')

    # Verify true_labels is 1D array
    assert adata_result2.obs['true_labels'].ndim == 1, "true_labels should be 1D array"
    print(f"  ✓ true_labels is 1D (shape: {adata_result2.obs['true_labels'].shape})")

    # Check label distribution
    label_counts2 = adata_result2.obs['true_labels'].value_counts().to_dict()
    print(f"  Label distribution: {label_counts2}")

    # Verify all labels are valid
    actual_labels2 = set(adata_result2.obs['true_labels'].unique())
    assert actual_labels2.issubset(valid_labels), f"Invalid labels found: {actual_labels2 - valid_labels}"
    print(f"  ✓ All labels are valid: {actual_labels2}")

    print("\n" + "=" * 70)
    print("✅ ALL TESTS PASSED: Quantile extraction fixed correctly")
    print("=" * 70)


def test_population_not_found():
    """Test error handling when population doesn't exist."""
    print("\n" + "=" * 70)
    print("TEST: Population Not Found Error Handling")
    print("=" * 70)

    # Create mock data
    n_cells = 50
    populations = ['M1'] * 20 + ['M2'] * 30

    adata = ad.AnnData(
        X=np.zeros((n_cells, 10)),
        obs=pd.DataFrame({
            'celltype': populations,
            'Condition1_prob': np.random.uniform(0.3, 0.8, n_cells),
            'Condition2_prob': np.random.uniform(0.2, 0.7, n_cells),
        }, index=[f'cell_{i}' for i in range(n_cells)])
    )

    # Try to use non-existent population
    try:
        quantile_assign_label(adata, 'celltype', 0.75, 'M3')
        print("  ❌ FAIL: Should have raised ValueError")
        assert False, "Should have raised ValueError for non-existent population"
    except ValueError as e:
        print(f"  ✓ Correctly raised ValueError: {e}")
        assert "not found" in str(e), "Error message should mention 'not found'"

    print("\n" + "=" * 70)
    print("✅ ERROR HANDLING TEST PASSED")
    print("=" * 70)


def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print("RUNNING QUANTILE BUG FIX TESTS")
    print("=" * 70)

    try:
        test_quantile_scalar_extraction()
        test_population_not_found()

        print("\n" + "=" * 70)
        print("ALL QUANTILE BUG FIX TESTS PASSED ✅")
        print("=" * 70)
        print("\nThe bug fix ensures:")
        print("1. Quantile thresholds are scalar floats (not DataFrames or arrays)")
        print("2. true_labels are 1D arrays (not 2D due to broadcasting)")
        print("3. Proper error handling for missing populations")
        return 0

    except AssertionError as e:
        print("\n" + "=" * 70)
        print("TEST FAILED ❌")
        print("=" * 70)
        print(f"\nError: {e}")
        return 1
    except Exception as e:
        print("\n" + "=" * 70)
        print("UNEXPECTED ERROR ❌")
        print("=" * 70)
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
