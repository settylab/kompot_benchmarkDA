"""
Test script to verify the fix for set_relevant_prob discontinuities.

This test demonstrates that the fixed implementation produces smooth gradients
without artificial discontinuities at cluster boundaries.
"""
import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add project root to path (two levels up from tests/unit/)
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from lib.condition_prob_centroid import create_enrichment_scores, set_relevant_prob
from lib.helper_functions import scale
import anndata as ad


def test_smooth_gradients():
    """
    Test that probabilities are smooth across cluster boundaries.
    """
    print("="*70)
    print("TESTING FIXED set_relevant_prob IMPLEMENTATION")
    print("="*70)

    # Create mock data: 3 populations (M1, M2, M3)
    # Enriching M1=0.95 and M2=0.75

    # Simulated fuzzy membership (sigmoid_fuzzy_weights after sigmoid transformation)
    # Each column = proximity to a centroid
    sigmoid_fuzzy_weights = pd.DataFrame({
        'M1': [0.95, 0.70, 0.40, 0.25, 0.20, 0.10, 0.05],
        'M2': [0.05, 0.15, 0.30, 0.20, 0.60, 0.80, 0.90],
        'M3': [0.10, 0.20, 0.50, 0.85, 0.40, 0.30, 0.15],
    }, index=[
        'cell_M1_center',
        'cell_M1_medium',
        'cell_M1_edge',
        'cell_M3_near_M1',  # M3 cell close to M1!
        'cell_M3_center',
        'cell_M3_near_M2',  # M3 cell close to M2!
        'cell_M2_edge',
    ])

    # Create mock AnnData
    adata = ad.AnnData(X=np.zeros((7, 10)))
    adata.obs_names = sigmoid_fuzzy_weights.index
    adata.obs['population'] = ['M1', 'M1', 'M1', 'M3', 'M3', 'M3', 'M2']

    print("\nInput sigmoid_fuzzy_weights (fuzzy membership to centroids):")
    print(sigmoid_fuzzy_weights.round(3))

    # Step 1: Create enrichment scores
    enr_scores = create_enrichment_scores(
        pop=['M1', 'M2'],
        pop_enr=[0.95, 0.75],
        populations=sigmoid_fuzzy_weights.columns
    )

    print("\nEnrichment scores:")
    print(enr_scores)

    # Step 2: Compute final probabilities (REFACTORED - single function)
    cond_prob = set_relevant_prob(
        sigmoid_fuzzy_weights=sigmoid_fuzzy_weights,
        enr_scores=enr_scores,
        pop=['M1', 'M2'],
        adata=adata,
        pop_column='population'
    )

    print("\nFinal cond_probability (FIXED):")
    for cell, prob in cond_prob.items():
        label = adata.obs.loc[cell, 'population']
        print(f"  {cell:20s} (label={label}): {prob:.3f}")

    # Check for smoothness
    print("\n" + "="*70)
    print("CHECKING FOR SMOOTH GRADIENTS:")
    print("="*70)

    # Check M1/M3 boundary
    m1_edge = cond_prob['cell_M1_edge']
    m3_near_m1 = cond_prob['cell_M3_near_M1']
    diff_1 = abs(m1_edge - m3_near_m1)

    print(f"\nM1/M3 boundary:")
    print(f"  cell_M1_edge (M1 label):     {m1_edge:.3f}")
    print(f"  cell_M3_near_M1 (M3 label):  {m3_near_m1:.3f}")
    print(f"  Difference: {diff_1:.3f}")

    # Check M3/M2 boundary
    m3_near_m2 = cond_prob['cell_M3_near_M2']
    m2_edge = cond_prob['cell_M2_edge']
    diff_2 = abs(m3_near_m2 - m2_edge)

    print(f"\nM3/M2 boundary:")
    print(f"  cell_M3_near_M2 (M3 label):  {m3_near_m2:.3f}")
    print(f"  cell_M2_edge (M2 label):     {m2_edge:.3f}")
    print(f"  Difference: {diff_2:.3f}")

    # Verify smooth gradient
    print("\n" + "="*70)
    print("GRADIENT VISUALIZATION:")
    print("="*70)
    print("\nM1 → M3 → M2:")
    for cell in sigmoid_fuzzy_weights.index:
        label = adata.obs.loc[cell, 'population']
        prob = cond_prob[cell]
        bar = '█' * int(prob * 50)
        print(f"  {cell:20s} ({label}): {prob:.3f} {bar}")

    # Success criteria
    print("\n" + "="*70)
    print("SUCCESS CRITERIA:")
    print("="*70)

    # All probabilities should be >= 0.5
    all_above_neutral = (cond_prob >= 0.5).all()
    print(f"✓ All probabilities ≥ 0.5: {all_above_neutral}")

    # Enriched populations should have min probability = 0.5
    m1_cells = cond_prob[adata.obs['population'] == 'M1']
    m2_cells = cond_prob[adata.obs['population'] == 'M2']
    enriched_min = min(m1_cells.min(), m2_cells.min())
    print(f"✓ Minimum in enriched populations: {enriched_min:.3f} (should be ~0.5)")

    # Check smoothness (small differences at boundaries)
    max_diff = max(diff_1, diff_2)
    print(f"✓ Maximum difference at boundaries: {max_diff:.3f} (should be smooth)")

    if all_above_neutral and abs(enriched_min - 0.5) < 0.01 and max_diff < 0.2:
        print("\n✅ TEST PASSED: Smooth gradients without discontinuities!")
    else:
        print("\n❌ TEST FAILED: Check the output above")

    return cond_prob


if __name__ == "__main__":
    test_smooth_gradients()
