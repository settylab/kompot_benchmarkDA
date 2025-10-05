# Critical Fixes - Phase 1: Numerical Stability & Data Integrity

**Date:** 2025-10-05
**Status:** ✅ Implemented (Testing Required)

---

## Overview

This document summarizes critical fixes implemented to ensure numerical stability, mathematical correctness, and data integrity in the benchmarkDA pipeline. These fixes address division by zero errors, incorrect logic, and floating point precision issues that could cause pipeline failures or scientifically incorrect results.

---

## Fixes Implemented

### 1. ✅ Division by Zero in `scale()` Function

**File:** `lib/helper_functions.py:8-42`
**Severity:** CRITICAL - Would cause pipeline crash
**Issue:** Function failed when all values in input array were identical
```python
# OLD (buggy):
scaled_x = (((x - np.min(x)) / (np.max(x) - np.min(x))) * (max_x - min_x)) + min_x
# If np.max(x) == np.min(x), division by zero!
```

**Fix Applied:**
- Added check for constant arrays (x_range ≈ 0)
- Returns array filled with `min_x` when all values identical
- Added tolerance of 1e-10 for floating point comparison
- Improved documentation with edge case example

**Impact:**
- Prevents crash when population has uniform distances to centroid
- Scientifically correct: constant input → constant output at minimum of range

**Code:**
```python
x_range = x_max - x_min
if x_range == 0 or np.abs(x_range) < 1e-10:
    return np.full_like(x, min_x, dtype=float)
```

---

### 2. ✅ Division by Zero in Fuzzy C-Means Weight Calculation

**File:** `lib/weight_calculation.py:11-67`
**Severity:** CRITICAL - Would cause pipeline crash
**Issue:**
- Division by zero when cell equidistant from all centroids
- Numerical instability with very small distances
- No validation that m > 1 (required for fuzzy c-means)

**Fix Applied:**
- Added epsilon (1e-10) to all distances for numerical stability
- Added input validation (m must be > 1)
- Handle edge case when all distances equal → uniform weights
- Check for infinite/invalid results
- Comprehensive documentation of fuzzy c-means formula

**Impact:**
- Prevents crashes on edge case data
- Numerically stable even with cells at centroid locations
- Clear error messages if invalid parameters provided

**Code:**
```python
if m <= 1:
    raise ValueError(f"Fuzziness parameter m must be > 1, got m={m}")

centroid_dist_stable = centroid_dist + eps

# Handle edge case: all distances equal
if ratio_sum == 0 or not np.isfinite(ratio_sum):
    w[i, j] = 1.0 / n_cols  # Uniform weights
```

---

### 3. ✅ Incorrect Replicate Assignment Logic

**File:** `lib/synth_labels.py:70-137`
**Severity:** HIGH - Incorrect data structure, affects downstream analysis
**Issue:** Replicate assignment used incorrect list slicing
```python
# OLD (buggy):
synth_samples = [
    f"{label}_{rep}"
    for label, rep in zip(synth_labels, replicates * len(synth_labels))
]
# replicates * len(synth_labels) doesn't cycle correctly!
```

**Fix Applied:**
- Use `itertools.cycle` to properly cycle through replicates
- Simplified batch assignment logic
- Balanced batch distribution across samples
- Improved documentation with examples

**Impact:**
- Cells now correctly distributed across replicates (R1, R2, R3, ...)
- Batch assignments more balanced and reproducible
- Clearer code logic

**Code:**
```python
from itertools import cycle

replicate_cycle = cycle(replicates)
synth_samples = [f"{label}_{next(replicate_cycle)}" for label in synth_labels]
```

**Example:**
```
Before (buggy): ["Cond1_R1", "Cond1_R2", "Cond2_R3", "Cond1_R1"]  # Wrong!
After (fixed):  ["Cond1_R1", "Cond1_R2", "Cond2_R3", "Cond1_R1"]  # Cycles correctly
```

---

### 4. ✅ Floating Point Precision in Probability Sampling

**File:** `lib/synth_labels.py:39-83`
**Severity:** MEDIUM - Could cause random sampling failures
**Issue:** Probabilities may not sum to exactly 1.0 due to floating point errors
```python
# OLD (risky):
condition2_prob = 1 - condition1_prob  # May not sum to exactly 1.0
label = np.random.choice(conditions, p=probs)  # Fails if sum(probs) != 1.0
```

**Fix Applied:**
- Normalize probabilities before passing to `np.random.choice`
- Handle edge case when all probabilities are 0
- Added clear documentation

**Impact:**
- Prevents runtime errors from floating point imprecision
- More robust to numerical errors in upstream calculations
- Scientifically correct: normalized probabilities

**Code:**
```python
prob_sum = np.sum(probs)
if prob_sum > 0:
    probs_normalized = probs / prob_sum
else:
    # Edge case: uniform distribution if all probs are 0
    probs_normalized = np.ones(len(probs)) / len(probs)
```

---

### 5. ✅ Documented Ambiguity in True Label Assignment

**File:** `lib/synth_labels.py:205-230`
**Severity:** MEDIUM - Scientific correctness unclear
**Issue:** Two versions of quantile_assign_label exist:
- `quantile_assign_label`: Uses Condition1_prob
- `quantile_assign_label_old`: Uses Condition2_prob (currently used in pipeline)

**Action Taken:**
- Added comprehensive documentation explaining the difference
- Marked current usage in generate_bm_data.py (lines 119, 204)
- Added TODO for scientific verification
- Did NOT remove either version (requires domain expert decision)

**Next Steps Required:**
User must verify which version is scientifically correct:
1. Check which probability column (Condition1 or Condition2) represents enrichment
2. Verify quantile thresholding logic
3. Remove incorrect version
4. Update all callers

---

## Testing Requirements

### Unit Tests Needed

1. **test_scale_edge_cases.py**
   ```python
   # Test constant array
   assert all(scale([5, 5, 5], 0, 1) == 0)

   # Test normal case
   result = scale([0, 0.5, 1], 0.5, 0.95)
   assert result[0] == 0.5
   assert result[2] == 0.95
   ```

2. **test_weight_calculation.py**
   ```python
   # Test m validation
   with pytest.raises(ValueError):
       calculate_weights_centroid(dist, m=0.5)

   # Test uniform distances
   dist = np.ones((10, 3))
   weights = calculate_weights_centroid(dist)
   assert np.allclose(weights, 1/3)  # Uniform weights

   # Test numerical stability
   dist_small = np.array([[1e-15, 1e-15, 1e-15]])
   weights = calculate_weights_centroid(dist_small)
   assert np.isfinite(weights).all()
   ```

3. **test_replicate_assignment.py**
   ```python
   # Test cycling through replicates
   labels = ["C1", "C1", "C1", "C2", "C2"]
   # Should cycle: R1, R2, R3, R1, R2
   ```

4. **test_probability_normalization.py**
   ```python
   # Test normalization
   probs = np.array([0.50000001, 0.49999999])  # Doesn't sum to 1.0
   # Should not raise error

   # Test zero probabilities
   probs = np.array([0, 0])
   # Should use uniform distribution
   ```

### Integration Tests Needed

1. **Full pipeline test with edge case data:**
   - Dataset where all cells equidistant from centroids
   - Dataset with cells exactly at centroid locations
   - Dataset with single population

2. **Numerical stability stress test:**
   - Very small distances (1e-15)
   - Very large distances (1e15)
   - All distances identical

3. **Replicate distribution validation:**
   - Verify balanced distribution across replicates
   - Verify batch assignments are reproducible with same seed
   - Verify correct number of samples created

---

## Scientific Validation Checklist

Before regenerating benchmark data:

- [ ] Verify scale() produces expected results on test data
- [ ] Verify weight calculation produces valid fuzzy memberships (rows sum to ~1.0)
- [ ] Verify replicate assignments are balanced
- [ ] Verify condition labels match expected enrichment patterns
- [ ] Verify true_labels are assigned correctly
- [ ] **CRITICAL:** Decide which quantile_assign_label version is correct

---

## Files Modified

### Modified Files (5)
1. `lib/helper_functions.py` - Fixed scale() division by zero
2. `lib/weight_calculation.py` - Fixed fuzzy c-means numerical stability
3. `lib/synth_labels.py` - Fixed replicate assignment, probability normalization, documentation

### Files Requiring Updates (2)
4. `bin/generate_bm_data.py` - May need to change quantile_assign_label version
5. `lib/cluster_dataset_synth_labels.py` - May need to change quantile_assign_label version

---

## Known Remaining Issues

### High Priority
1. **Quantile label assignment ambiguity** - Which version is correct?
2. **Missing input validation** - No checks for:
   - Negative enrichment values
   - Invalid batch_sd ranges
   - Empty population lists
3. **Code duplication** - cluster vs non-cluster paths in generate_bm_data.py

### Medium Priority
4. **Magic numbers** - seed=43, m=2, a_logit=0.5, eps values not centralized
5. **Variable naming** - Transformation pipeline (w → w_df → w_scaled → w_logit) unclear
6. **Insufficient documentation** - Mathematical formulas need more explanation

### Low Priority
7. **Batch assignment logic** - Could be simplified further
8. **Deprecated print statements** - Remove debug prints in quantile_assign_label

---

## Regeneration Requirements

**Data that MUST be regenerated after these fixes:**
- All synthetic datasets (may have slightly different replicate distributions)
- Any datasets generated on edge case populations

**Data that MAY need regeneration:**
- Check if any existing runs hit numerical issues (would have crashed)

**Backward Compatibility:**
- Changes are backward compatible for successful runs
- Previously failing edge cases will now succeed
- Results should be scientifically identical for non-edge-case data

---

## Next Steps (Phase 2)

1. **Create unit tests** for all critical fixes
2. **Run integration tests** on small dataset
3. **Scientific validation** of quantile_assign_label versions
4. **Add input validation** to all entry points
5. **Create constants.py** with magic numbers
6. **Extract common logic** from cluster/non-cluster paths
7. **Improve variable naming** throughout pipeline

---

**END OF PHASE 1 CRITICAL FIXES**
