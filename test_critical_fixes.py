"""
Test script for Phase 1 critical fixes.

Tests numerical stability improvements and bug fixes:
1. scale() function with constant arrays
2. calculate_weights_centroid() with edge cases
3. Replicate assignment cycling
4. Probability normalization
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from lib.helper_functions import scale
from lib.weight_calculation import calculate_weights_centroid
from lib.synth_labels import label_condition_and_rep_labels, label_condition_and_rep_other
import anndata as ad


def test_scale_constant_array():
    """Test that scale() handles constant arrays without division by zero."""
    print("=" * 70)
    print("TEST 1: scale() with constant array")
    print("=" * 70)

    # Test constant array
    x_constant = np.array([5.0, 5.0, 5.0, 5.0])
    result = scale(x_constant, 0.5, 0.95)

    print(f"Input (constant):  {x_constant}")
    print(f"Output:            {result}")
    print(f"Expected:          [0.5, 0.5, 0.5, 0.5]")

    # Should return array of min_x
    assert np.allclose(result, 0.5), f"Expected all 0.5, got {result}"
    print("✅ PASS: Constant array handled correctly\n")


def test_scale_normal():
    """Test that scale() still works correctly on normal arrays."""
    print("=" * 70)
    print("TEST 2: scale() with normal array")
    print("=" * 70)

    x_normal = np.array([0.0, 0.5, 1.0])
    result = scale(x_normal, 0.5, 0.95)

    print(f"Input:    {x_normal}")
    print(f"Output:   {result}")
    print(f"Expected: [0.5, 0.725, 0.95]")

    assert np.isclose(result[0], 0.5), f"Expected 0.5, got {result[0]}"
    assert np.isclose(result[1], 0.725), f"Expected 0.725, got {result[1]}"
    assert np.isclose(result[2], 0.95), f"Expected 0.95, got {result[2]}"
    print("✅ PASS: Normal scaling works correctly\n")


def test_weight_calculation_uniform_distances():
    """Test fuzzy c-means with uniform distances."""
    print("=" * 70)
    print("TEST 3: calculate_weights_centroid() with uniform distances")
    print("=" * 70)

    # All distances equal
    dist = np.ones((5, 3))
    weights = calculate_weights_centroid(dist, m=2)

    print(f"Input distances (all 1.0):")
    print(dist)
    print(f"\nOutput weights:")
    print(weights)
    print(f"Expected: All weights ≈ 1/3 = 0.333...")

    # Should give uniform weights (each cell belongs equally to all centroids)
    expected = 1.0 / 3
    assert np.allclose(weights, expected, rtol=0.01), f"Expected uniform weights ~{expected}, got {weights}"

    # Check rows sum to ~1.0
    row_sums = weights.sum(axis=1)
    print(f"\nRow sums: {row_sums}")
    assert np.allclose(row_sums, 1.0, rtol=0.01), f"Expected row sums ~1.0, got {row_sums}"
    print("✅ PASS: Uniform distances handled correctly\n")


def test_weight_calculation_numerical_stability():
    """Test fuzzy c-means with very small distances."""
    print("=" * 70)
    print("TEST 4: calculate_weights_centroid() with tiny distances")
    print("=" * 70)

    # Very small distances (near zero)
    dist = np.array([[1e-10, 2e-10, 3e-10],
                     [1e-15, 1e-15, 1e-15]])

    weights = calculate_weights_centroid(dist, m=2)

    print(f"Input distances (very small):")
    print(dist)
    print(f"\nOutput weights:")
    print(weights)

    # Should not have NaN or Inf
    assert np.isfinite(weights).all(), f"Weights contain NaN or Inf: {weights}"

    # Rows should sum to ~1.0
    row_sums = weights.sum(axis=1)
    print(f"Row sums: {row_sums}")
    assert np.allclose(row_sums, 1.0, rtol=0.1), f"Expected row sums ~1.0, got {row_sums}"
    print("✅ PASS: Numerical stability maintained\n")


def test_weight_calculation_m_validation():
    """Test that m parameter is validated."""
    print("=" * 70)
    print("TEST 5: calculate_weights_centroid() m parameter validation")
    print("=" * 70)

    dist = np.ones((5, 3))

    # m must be > 1
    try:
        weights = calculate_weights_centroid(dist, m=0.5)
        assert False, "Should have raised ValueError for m <= 1"
    except ValueError as e:
        print(f"✅ Correctly raised ValueError: {e}")

    try:
        weights = calculate_weights_centroid(dist, m=1.0)
        assert False, "Should have raised ValueError for m <= 1"
    except ValueError as e:
        print(f"✅ Correctly raised ValueError: {e}")

    # m > 1 should work
    weights = calculate_weights_centroid(dist, m=2.0)
    print(f"✅ m=2.0 works correctly")
    print("✅ PASS: m parameter validation works\n")


def test_replicate_cycling():
    """Test that replicate assignment cycles correctly."""
    print("=" * 70)
    print("TEST 6: Replicate assignment cycling")
    print("=" * 70)

    # Create mock adata
    n_cells = 7
    adata = ad.AnnData(X=np.zeros((n_cells, 10)))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]

    # Assign synthetic labels (alternating conditions)
    synth_labels = ["Condition1", "Condition2", "Condition1", "Condition2",
                    "Condition1", "Condition2", "Condition1"]
    adata.obs["synth_labels"] = synth_labels

    # Test with 3 replicates
    adata = label_condition_and_rep_other(adata, n_replicates=3, n_batches=2, seed=42)

    print("Synth labels:  ", synth_labels)
    print("Synth samples: ", list(adata.obs["synth_samples"]))

    # Check that replicates cycle
    samples = list(adata.obs["synth_samples"])

    # Should cycle through R1, R2, R3
    assert "R1" in samples[0], f"First cell should be R1, got {samples[0]}"
    assert "R2" in samples[1], f"Second cell should be R2, got {samples[1]}"
    assert "R3" in samples[2], f"Third cell should be R3, got {samples[2]}"

    # Check all cells have synth_samples and synth_batches
    assert all(pd.notna(adata.obs["synth_samples"])), "Some cells missing synth_samples"
    assert all(pd.notna(adata.obs["synth_batches"])), "Some cells missing synth_batches"

    print("✅ PASS: Replicate cycling works correctly\n")


def test_probability_normalization():
    """Test that probabilities are normalized before sampling."""
    print("=" * 70)
    print("TEST 7: Probability normalization")
    print("=" * 70)

    # Create mock adata
    n_cells = 100
    adata = ad.AnnData(X=np.zeros((n_cells, 10)))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]

    # Create probabilities that don't sum to exactly 1.0 (floating point errors)
    cond_prob = pd.DataFrame({
        "Condition1": [0.50000001] * n_cells,
        "Condition2": [0.49999999] * n_cells
    })

    print(f"Probability sum for first cell: {cond_prob.iloc[0].sum()}")
    print(f"(Should not be exactly 1.0 due to floating point)")

    # This should not raise an error
    try:
        adata = label_condition_and_rep_labels(adata, cond_prob, seed=42)
        print("✅ No error raised for imperfect probability sum")
    except Exception as e:
        assert False, f"Should not raise error, got: {e}"

    # Check that labels were assigned
    assert "synth_labels" in adata.obs.columns, "synth_labels not created"
    assert all(pd.notna(adata.obs["synth_labels"])), "Some cells missing labels"

    print("✅ PASS: Probability normalization works\n")


def test_zero_probabilities():
    """Test edge case where all probabilities are zero."""
    print("=" * 70)
    print("TEST 8: Zero probabilities edge case")
    print("=" * 70)

    # Create mock adata
    n_cells = 10
    adata = ad.AnnData(X=np.zeros((n_cells, 10)))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]

    # All probabilities are zero (edge case)
    cond_prob = pd.DataFrame({
        "Condition1": [0.0] * n_cells,
        "Condition2": [0.0] * n_cells
    })

    print(f"All probabilities are zero")

    # Should fall back to uniform distribution
    try:
        adata = label_condition_and_rep_labels(adata, cond_prob, seed=42)
        print("✅ No error raised for zero probabilities")
    except Exception as e:
        assert False, f"Should not raise error, got: {e}"

    # Should have assigned labels (uniformly)
    assert "synth_labels" in adata.obs.columns, "synth_labels not created"
    assert all(pd.notna(adata.obs["synth_labels"])), "Some cells missing labels"

    print(f"Labels assigned: {adata.obs['synth_labels'].value_counts().to_dict()}")
    print("✅ PASS: Zero probabilities handled correctly\n")


def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print("TESTING CRITICAL FIXES - PHASE 1")
    print("=" * 70 + "\n")

    try:
        test_scale_constant_array()
        test_scale_normal()
        test_weight_calculation_uniform_distances()
        test_weight_calculation_numerical_stability()
        test_weight_calculation_m_validation()
        test_replicate_cycling()
        test_probability_normalization()
        test_zero_probabilities()

        print("\n" + "=" * 70)
        print("ALL TESTS PASSED ✅")
        print("=" * 70)
        print("\nCritical fixes are working correctly!")
        print("Next steps:")
        print("1. Run full pipeline test on small dataset")
        print("2. Decide which quantile_assign_label version to use")
        print("3. Regenerate benchmark data if needed")
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
