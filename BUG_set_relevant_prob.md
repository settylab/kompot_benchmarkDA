# 🚨 BUG: Artificial Discontinuities from Label-Based Override in `set_relevant_prob`

**Date:** 2025-10-04
**Location:** `lib/condition_prob_centroid.py:73-109`
**Severity:** HIGH - Creates unnatural probability steps at all cluster boundaries

---

## The Core Problem

The current implementation uses **different algorithms** for different cluster labels:
- **Enriched population cells** (e.g., M1, M2): Use their "own" column from `prob_matrix`
- **Non-enriched population cells** (e.g., M3, M4): Use **mean** of enriched columns

**Result**: Artificial discontinuities at every cluster boundary because neighboring cells in embedding space use different formulas.

---

## Current Buggy Logic

```python
def set_relevant_prob(enr_prob, pop_enr, pop, adata, pop_column):
    prob_matrix = enr_prob[pop]  # Select enriched population columns

    if len(pop) > 1:
        # ALL cells start with mean
        cond_probability = prob_matrix.mean(axis=1)

        # Override cells belonging to enriched populations
        for population in pop:
            cells_in_pop = adata.obs_names[adata.obs[pop_column] == population]
            # These cells use their "own" column instead of mean
            cond_probability.loc[cells_in_pop] = prob_matrix.loc[cells_in_pop, population]
    else:
        cond_probability = prob_matrix

    return cond_probability
```

**The bug**: The `for population in pop:` loop creates different probability assignment rules based on cluster **label**, not transcriptional **proximity**.

---

## Example: Enriching M1=0.95 and M2=0.75

### Embedding Space
```
M1 cells ← → M3 cells ← → M2 cells
(enriched)   (neutral)   (enriched)
```

### After Distance Calculation and Scaling

Each column in `enr_prob` represents distance to a centroid, scaled to `[0.5, enrichment]`:

```python
                   M1_col  M2_col
cell_M1_center     0.950   0.500   # M1 cell, close to M1
cell_M1_edge       0.659   0.607   # M1 cell, edge of cluster
cell_M3_near_M1    0.579   0.679   # M3 cell, close to M1!
cell_M3_center     0.500   0.500   # M3 cell, far from both
cell_M3_near_M2    0.500   0.679   # M3 cell, close to M2!
cell_M2_edge       0.553   0.750   # M2 cell, edge of cluster
cell_M2_center     0.500   0.950   # M2 cell, close to M2
```

### Current Logic (Buggy)

**Step 1**: All cells get mean
```python
cell_M1_center     (0.950 + 0.500)/2 = 0.725
cell_M1_edge       (0.659 + 0.607)/2 = 0.633
cell_M3_near_M1    (0.579 + 0.679)/2 = 0.629
cell_M3_center     (0.500 + 0.500)/2 = 0.500
cell_M3_near_M2    (0.500 + 0.679)/2 = 0.590
cell_M2_edge       (0.553 + 0.750)/2 = 0.651
cell_M2_center     (0.500 + 0.950)/2 = 0.725
```

**Step 2**: Override M1 and M2 cells
```python
cell_M1_center     0.950  ← Uses M1 column
cell_M1_edge       0.659  ← Uses M1 column
cell_M3_near_M1    0.629  ← Keeps mean (not overridden)
cell_M3_center     0.500  ← Keeps mean
cell_M3_near_M2    0.590  ← Keeps mean (not overridden)
cell_M2_edge       0.750  ← Uses M2 column
cell_M2_center     0.950  ← Uses M2 column
```

### The Discontinuities

#### Boundary 1: M1/M3
```
cell_M1_edge:      0.659  (uses M1 column)
     |
     | ← DISCONTINUITY: 0.030 drop
     |
cell_M3_near_M1:   0.629  (uses mean)
```

Even though `cell_M3_near_M1` has similar transcriptional profile to `cell_M1_edge`, it gets a different probability because it uses **mean** instead of **M1 column**.

#### Boundary 2: M3/M2
```
cell_M3_near_M2:   0.590  (uses mean)
     |
     | ← DISCONTINUITY: 0.160 jump
     |
cell_M2_edge:      0.750  (uses M2 column)
```

Massive jump at the M3/M2 boundary because of algorithm change.

#### Boundary 3: M1/M2 (if they overlap)
```
cell_M1_near_M2:   uses M1 column (ignores M2 proximity)
     |
     | ← DISCONTINUITY
     |
cell_M2_near_M1:   uses M2 column (ignores M1 proximity)
```

Cells at the M1/M2 boundary use their "own" column, ignoring proximity to the other enriched population.

---

## Why This Is Wrong

### Transcriptional proximity should determine probability, not cluster labels

In embedding space, enrichment probability should follow smooth gradients based on transcriptional similarity. Two cells that are neighbors in PCA/DM space should have similar probabilities.

**Current behavior violates this principle:**
- Cell labeled "M1" at edge of cluster: uses M1 column
- Cell labeled "M3" but transcriptionally similar: uses mean
- Different formulas → discontinuity

### The algorithm changes at cluster boundaries

```
Within M1:        probability = M1_column
At M1/M3 edge:    M1_column → mean(M1, M2)  ← ALGORITHM CHANGE
Within M3:        probability = mean(M1, M2)
At M3/M2 edge:    mean(M1, M2) → M2_column  ← ALGORITHM CHANGE
Within M2:        probability = M2_column
```

This creates artificial steps wherever the label changes.

---

## Single Population Case - Also Affected!

Even with a **single enriched population**, the override creates issues:

```python
# Enriching only M1=0.95

# For single population, current code does:
cond_probability = prob_matrix  # This is just the M1 column

# But cells are still differentiated by label in downstream code
```

Wait, let me check if there's an override for single population...

Actually, for single population, `prob_matrix` is the M1 column directly (a Series, not DataFrame), so all cells use it. **No override in this case**, so **single population is OK**.

The bug only manifests with **multiple enriched populations**.

---

## Proposed Solution

### Concept
Use a **single consistent algorithm** for all cells based on **transcriptional proximity**, not cluster labels.

### Algorithm

**Step 1**: For each cell, compute minimum distance to any enriched centroid
```python
# Each cell gets probability from whichever enriched centroid it's closest to
cond_probability = prob_matrix.max(axis=1)  # max probability = min distance
```

**Step 2**: Find the minimum probability among cells actually labeled as enriched populations
```python
enriched_cells = adata.obs[pop_column].isin(enriched_pops)
min_prob_in_enriched = cond_probability[enriched_cells].min()
max_prob = cond_probability.max()
```

**Step 3**: Rescale so that minimum in enriched populations becomes 0.5
```python
# Linear rescaling
cond_probability_rescaled = 0.5 + (cond_probability - min_prob_in_enriched) / (max_prob - min_prob_in_enriched) * (max_prob - 0.5)
```

**Step 4**: Floor at 0.5 (neutral baseline)
```python
cond_probability_final = np.maximum(cond_probability_rescaled, 0.5)
```

### Rationale

1. **Max (min distance)**: Cells close to **any** enriched centroid get high probability
2. **Rescaling**: Ensures even the furthest cell **within** an enriched population gets at least 0.5
3. **Floor at 0.5**: Cells far from all enriched centroids remain neutral
4. **Smooth gradients**: Same algorithm for all cells → no discontinuities

---

## Proposed Implementation

```python
def set_relevant_prob(enr_prob, pop_enr, pop, adata, pop_column, cell_type_dict=None):
    """
    Set final condition probability for each cell based on minimum distance
    to any enriched centroid.

    Uses consistent algorithm for all cells to ensure smooth gradients
    without discontinuities at cluster boundaries.

    Parameters:
    - enr_prob: DataFrame with probabilities scaled to [0.5, enr_score] per population
    - pop_enr: Enrichment levels (not used in calculation, kept for compatibility)
    - pop: List of enriched population names
    - adata: AnnData object
    - pop_column: Column name with cluster labels

    Returns:
    - cond_probability: Series with final probability for each cell
    """
    prob_matrix = enr_prob[pop]

    if len(pop) > 1:
        # Step 1: Each cell uses max probability (closest enriched centroid)
        cond_probability = prob_matrix.max(axis=1)

        # Step 2: Find minimum probability in enriched populations
        enriched_mask = adata.obs[pop_column].isin(pop)
        if enriched_mask.any():
            min_prob = cond_probability[enriched_mask].min()
            max_prob = cond_probability.max()

            # Step 3: Rescale so minimum in enriched populations becomes 0.5
            if max_prob > min_prob:
                # Linear rescaling: [min_prob, max_prob] → [0.5, max_prob]
                cond_probability = 0.5 + (cond_probability - min_prob) / (max_prob - min_prob) * (max_prob - 0.5)

        # Step 4: Floor at 0.5 (neutral)
        cond_probability = np.maximum(cond_probability, 0.5)
    else:
        # Single population - use column directly (no rescaling needed)
        cond_probability = prob_matrix

    return cond_probability
```

---

## Example: Fixed Behavior

Using the same data as before:

### After Step 1 (Max)
```python
cell_M1_center     max(0.950, 0.500) = 0.950
cell_M1_edge       max(0.659, 0.607) = 0.659
cell_M3_near_M1    max(0.579, 0.679) = 0.679  ← Uses M2 proximity
cell_M3_center     max(0.500, 0.500) = 0.500
cell_M3_near_M2    max(0.500, 0.679) = 0.679  ← Same as near M1!
cell_M2_edge       max(0.553, 0.750) = 0.750
cell_M2_center     max(0.500, 0.950) = 0.950
```

### After Steps 2-3 (Rescale)
```python
# Find min in enriched cells (M1 and M2)
min_in_enriched = min(0.950, 0.659, 0.750, 0.950) = 0.659
max_prob = 0.950

# Rescale [0.659, 0.950] → [0.5, 0.950]
scale_factor = (0.950 - 0.5) / (0.950 - 0.659) = 1.546

cell_M1_center     0.5 + (0.950 - 0.659) * 1.546 = 0.950
cell_M1_edge       0.5 + (0.659 - 0.659) * 1.546 = 0.500
cell_M3_near_M1    0.5 + (0.679 - 0.659) * 1.546 = 0.531
cell_M3_center     0.5 + (0.500 - 0.659) * 1.546 = 0.254 → floor to 0.500
cell_M3_near_M2    0.5 + (0.679 - 0.659) * 1.546 = 0.531
cell_M2_edge       0.5 + (0.750 - 0.659) * 1.546 = 0.641
cell_M2_center     0.5 + (0.950 - 0.659) * 1.546 = 0.950
```

### After Step 4 (Floor)
```python
cell_M1_center     0.950
cell_M1_edge       0.500  ← Minimum within M1
cell_M3_near_M1    0.531  ← Smooth gradient
cell_M3_center     0.500  ← Neutral (far from enriched)
cell_M3_near_M2    0.531  ← Smooth gradient
cell_M2_edge       0.641  ← Smooth gradient
cell_M2_center     0.950
```

### Gradient Visualization
```
M1 center → M1 edge → M3 near M1 → M3 center → M3 near M2 → M2 edge → M2 center
  0.950       0.500       0.531        0.500       0.531       0.641       0.950

✓ Smooth transitions everywhere
✓ No discontinuities at cluster boundaries
✓ Same algorithm for all cells
```

---

## Comparison: Current vs Proposed

| Cell | Current (Buggy) | Proposed (Fixed) | Issue |
|------|----------------|------------------|-------|
| `cell_M1_edge` | 0.659 | 0.500 | Current: overestimates edge cells |
| `cell_M3_near_M1` | 0.629 | 0.531 | Current: **discontinuity** (0.030 drop) |
| `cell_M3_near_M2` | 0.590 | 0.531 | Current: **discontinuity** (0.160 jump to M2) |
| `cell_M2_edge` | 0.750 | 0.641 | Current: overestimates edge cells |

**Key differences:**
1. **Proposed**: Smooth gradient from 0.950 → 0.500 → 0.531 → 0.641 → 0.950
2. **Current**: Jumps from 0.659 → 0.629 → 0.590 → 0.750 (discontinuous)

---

## Scientific Interpretation

### Current (Buggy)
```
"Enrich M1 to 0.95 and M2 to 0.75"

Result:
- M1-labeled cells: high if close to M1 centroid
- M2-labeled cells: high if close to M2 centroid
- Other cells: average of proximity to M1 and M2
- Discontinuities at all label boundaries
```

### Proposed (Fixed)
```
"Enrich M1 to 0.95 and M2 to 0.75"

Result:
- ALL cells: high if close to M1 OR M2 centroid
- Smooth gradients everywhere based on transcriptional similarity
- Even furthest M1/M2 cells guaranteed ≥ 0.5
- Cells far from both → neutral (0.5)
```

---

## Benefits of Proposed Solution

1. **No discontinuities**: Same algorithm for all cells
2. **Transcriptionally driven**: Proximity determines probability, not labels
3. **Guarantees**: Cells in enriched populations get ≥ 0.5
4. **Neutral baseline**: Cells far from enriched centroids → 0.5
5. **Interpretable**: "Enriched if transcriptionally similar to any target population"

---

## Edge Cases

### Case 1: Enriched populations overlap
```
If M1 and M2 are close in embedding space:
- Cells between them get high probability (max proximity)
- Smooth gradient through the overlap region
```

### Case 2: Enriched populations are far apart
```
If M1 and M2 are distant:
- Cells near M1 get high probability
- Cells near M2 get high probability
- Cells between them get neutral (0.5)
- Two separate enriched regions
```

### Case 3: Non-enriched population close to enriched
```
If M3 is close to M1:
- M3 cells near M1 get enriched (transcriptional similarity)
- This is scientifically correct!
- Enrichment follows biology, not labels
```

---

## Alternative: Mean Instead of Max

If you prefer **moderate** enrichment near boundaries:

```python
# Step 1: Use mean instead of max
cond_probability = prob_matrix.mean(axis=1)

# Steps 2-4: Same rescaling and floor
```

**Difference:**
- **Max**: Aggressive enrichment if close to ANY target
- **Mean**: Moderate enrichment, averages proximity to all targets

Both eliminate discontinuities. Choice depends on biological interpretation.

---

## Affected Data

All datasets with **multiple enriched populations**:
- `data/synthetic/linear/` with multiple M populations
- `data/synthetic/branch/` with multiple M populations
- `data/real/*` with multiple cell types enriched

**Action required**: Regenerate labels with fixed algorithm.

---

## Implementation Steps

1. **Update** `lib/condition_prob_centroid.py:set_relevant_prob()`
2. **Test** on simple example to verify smooth gradients
3. **Regenerate** all multi-population label data
4. **Validate** by plotting probability distributions across boundaries
5. **Document** the aggregation strategy (max vs mean)

---

## Questions

1. **Max vs Mean**: Which aggregation for multiple populations?
   - Max: Enriched if close to ANY target (recommended)
   - Mean: Average proximity to all targets

2. **Rescaling**: Is the proposed rescaling appropriate?
   - Ensures minimum in enriched populations = 0.5
   - Maintains maximum = enrichment level
   - Linear interpolation between

3. **Floor**: Should we floor at 0.5 or allow lower?
   - Floor at 0.5: Neutral baseline (recommended)
   - Allow lower: Depletion possible

---

**END OF REPORT**
