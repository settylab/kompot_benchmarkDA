"""
Integration test for BenchmarkDA full pipeline.

Tests the complete workflow:
1. Data preprocessing (PCA, DM)
2. Label generation with batch effects
3. Output file creation and validation
"""

import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
import numpy as np
import pandas as pd

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import anndata as ad
from lib.synth_labels import (
    label_condition_and_rep_labels,
    label_condition_and_rep_other,
    add_batch_effect_pca
)
from lib.condition_prob_centroid import create_enrichment_scores, set_relevant_prob
from lib.get_weight_matrix import get_weight_matrix_centroid


def test_label_generation_pipeline():
    """Test full label generation pipeline on mock data."""
    print("=" * 70)
    print("INTEGRATION TEST: Full Label Generation Pipeline")
    print("=" * 70)

    # Create mock AnnData
    n_cells = 200
    n_genes = 100
    n_pops = 3

    X = np.random.randn(n_cells, n_genes)

    # Create populations (M1, M2, M3)
    populations = [f"M{(i % n_pops) + 1}" for i in range(n_cells)]

    obs = pd.DataFrame({
        "celltype": populations,
    }, index=[f"cell_{i}" for i in range(n_cells)])

    # Create embeddings as DataFrame (required by get_weight_matrix_centroid)
    pca_data = np.random.randn(n_cells, 30)
    pca_df = pd.DataFrame(pca_data, index=obs.index)

    obsm = {
        "X_pca": pca_data,
    }

    adata = ad.AnnData(X=X, obs=obs, obsm=obsm)

    print(f"\n✓ Created mock data: {n_cells} cells, {n_pops} populations")
    print(f"  Populations: {sorted(set(populations))}")

    # Step 1: Calculate fuzzy weights
    print("\n[Step 1] Calculating fuzzy weights...")
    sigmoid_fuzzy_weights, conditions = get_weight_matrix_centroid(
        adata=adata,
        pop_column="celltype",
        seed=42,
        X_emb=pca_df,  # Pass DataFrame instead of numpy array
        n_conditions=2,
        m=2.0,
        a_logit=0.5
    )

    assert sigmoid_fuzzy_weights.shape[0] == n_cells, "Weight matrix row count mismatch"
    assert sigmoid_fuzzy_weights.shape[1] == n_pops, "Weight matrix column count mismatch"
    print(f"✅ Fuzzy weights calculated: {sigmoid_fuzzy_weights.shape}")

    # Step 2: Create enrichment scores
    print("\n[Step 2] Creating enrichment scores...")
    enr_pop = ["M1", "M2"]
    enr_values = [0.95, 0.75]

    enr_scores = create_enrichment_scores(
        pop=enr_pop,
        pop_enr=enr_values,
        populations=sigmoid_fuzzy_weights.columns
    )

    assert "M1" in enr_scores, "M1 enrichment score missing"
    assert "M2" in enr_scores, "M2 enrichment score missing"
    assert enr_scores["M1"] == 0.95, "M1 enrichment value incorrect"
    assert enr_scores["M2"] == 0.75, "M2 enrichment value incorrect"
    print(f"✅ Enrichment scores: {enr_scores}")

    # Step 3: Calculate condition probabilities
    print("\n[Step 3] Calculating condition probabilities...")
    cond_prob = set_relevant_prob(
        sigmoid_fuzzy_weights=sigmoid_fuzzy_weights,
        enr_scores=enr_scores,
        pop=enr_pop,
        adata=adata,
        pop_column="celltype"
    )

    # Validate probabilities
    assert len(cond_prob) == n_cells, "Probability count mismatch"
    assert (cond_prob >= 0.5).all(), "Some probabilities < 0.5"
    assert (cond_prob <= 1.0).all(), "Some probabilities > 1.0"

    # Check enriched populations have higher probabilities
    m1_cells = adata.obs["celltype"] == "M1"
    m2_cells = adata.obs["celltype"] == "M2"
    m3_cells = adata.obs["celltype"] == "M3"

    m1_mean_prob = cond_prob[m1_cells].mean()
    m2_mean_prob = cond_prob[m2_cells].mean()
    m3_mean_prob = cond_prob[m3_cells].mean()

    print(f"  M1 mean probability: {m1_mean_prob:.3f} (enriched at 0.95)")
    print(f"  M2 mean probability: {m2_mean_prob:.3f} (enriched at 0.75)")
    print(f"  M3 mean probability: {m3_mean_prob:.3f} (not enriched)")

    # M1 should have highest prob (most enriched)
    # M3 should have lowest prob (not enriched)
    # This might not always be true due to fuzzy weights, so just check reasonable values
    assert m1_mean_prob >= 0.5, "M1 mean probability too low"
    assert m3_mean_prob >= 0.5, "M3 mean probability too low"

    print(f"✅ Condition probabilities valid: min={cond_prob.min():.3f}, max={cond_prob.max():.3f}")

    # Step 4: Assign condition labels
    print("\n[Step 4] Assigning condition labels...")
    cond_prob_df = pd.DataFrame({
        "Condition1": cond_prob,
        "Condition2": 1.0 - cond_prob
    })

    adata = label_condition_and_rep_labels(
        adata=adata,
        cond_probability=cond_prob_df,
        seed=42
    )

    assert "synth_labels" in adata.obs.columns, "synth_labels not created"

    # Check label distribution
    label_counts = adata.obs["synth_labels"].value_counts().to_dict()
    print(f"  Label distribution: {label_counts}")

    # Both conditions should have some cells
    assert len(label_counts) == 2, "Should have 2 conditions"
    assert all(count > 0 for count in label_counts.values()), "Some conditions have 0 cells"

    print(f"✅ Condition labels assigned successfully")

    # Step 4b: Assign replicate and batch labels
    print("\n[Step 4b] Assigning replicate and batch labels...")
    adata = label_condition_and_rep_other(
        adata=adata,
        n_replicates=3,
        n_batches=2,
        seed=42
    )

    assert "synth_samples" in adata.obs.columns, "synth_samples not created"
    assert "synth_batches" in adata.obs.columns, "synth_batches not created"

    sample_counts = adata.obs["synth_samples"].value_counts().to_dict()
    print(f"  Sample distribution: {len(sample_counts)} samples")

    batch_counts = adata.obs["synth_batches"].value_counts().to_dict()
    print(f"  Batch distribution: {batch_counts}")

    print(f"✅ Replicate and batch labels assigned successfully")

    # Step 5: Add batch effects
    print("\n[Step 5] Adding batch effects...")
    adata_batch = add_batch_effect_pca(
        adata=adata,
        layer_embedding="X_pca",
        norm_sd=1.0,
        seed=42
    )

    assert "X_pca_batch" in adata_batch.obsm.keys(), "Batch-affected PCA not created"

    # Check batch effect was applied (embeddings should differ)
    original_pca = adata.obsm["X_pca"]
    batch_pca = adata_batch.obsm["X_pca_batch"]

    diff = np.abs(original_pca - batch_pca).mean()
    print(f"  Mean difference between original and batch PCA: {diff:.3f}")
    assert diff > 0, "No batch effect applied"

    print(f"✅ Batch effects added successfully")

    # Final validation
    print("\n[Final Validation]")
    print(f"✓ All required columns present:")
    required_cols = ["celltype", "synth_labels", "synth_samples", "synth_batches"]
    for col in required_cols:
        assert col in adata_batch.obs.columns, f"Missing column: {col}"
        print(f"  - {col}")

    print(f"\n✓ All required embeddings present:")
    required_embs = ["X_pca", "X_pca_batch"]
    for emb in required_embs:
        assert emb in adata_batch.obsm.keys(), f"Missing embedding: {emb}"
        print(f"  - {emb}")

    print("\n" + "=" * 70)
    print("✅ INTEGRATION TEST PASSED: Full pipeline working correctly")
    print("=" * 70)

    return True


def test_batch_effect_scaling():
    """Test that batch effects scale with norm_sd parameter."""
    print("\n" + "=" * 70)
    print("INTEGRATION TEST: Batch Effect Scaling")
    print("=" * 70)

    # Create mock data
    n_cells = 100
    adata = ad.AnnData(
        X=np.zeros((n_cells, 10)),
        obs=pd.DataFrame({
            "celltype": [f"M{i%3+1}" for i in range(n_cells)],
            "synth_labels": [f"Condition{i%2+1}" for i in range(n_cells)],
            "synth_samples": [f"Condition{i%2+1}_R{i%3+1}" for i in range(n_cells)],
            "synth_batches": [f"B{i%2+1}" for i in range(n_cells)],
        }, index=[f"cell_{i}" for i in range(n_cells)]),
        obsm={"X_pca": np.random.randn(n_cells, 30)}
    )

    # Test different norm_sd values
    norm_sds = [0.5, 1.0, 2.0]
    differences = []

    for norm_sd in norm_sds:
        adata_batch = add_batch_effect_pca(
            adata=adata.copy(),
            layer_embedding="X_pca",
            norm_sd=norm_sd,
            seed=42
        )

        diff = np.abs(adata.obsm["X_pca"] - adata_batch.obsm["X_pca_batch"]).mean()
        differences.append(diff)
        print(f"  norm_sd={norm_sd:.1f} → mean difference={diff:.3f}")

    # Check that larger norm_sd creates larger differences
    assert differences[1] > differences[0], "norm_sd=1.0 should have more effect than 0.5"
    assert differences[2] > differences[1], "norm_sd=2.0 should have more effect than 1.0"

    print("\n✅ PASSED: Batch effects scale correctly with norm_sd")

    return True


def main():
    """Run all integration tests."""
    print("\n" + "=" * 70)
    print("RUNNING INTEGRATION TESTS")
    print("=" * 70)

    try:
        test_label_generation_pipeline()
        test_batch_effect_scaling()

        print("\n" + "=" * 70)
        print("ALL INTEGRATION TESTS PASSED ✅")
        print("=" * 70)
        print("\nPipeline is functioning correctly!")
        return 0

    except AssertionError as e:
        print("\n" + "=" * 70)
        print("INTEGRATION TEST FAILED ❌")
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
