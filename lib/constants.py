"""
Constants for the benchmarkDA pipeline.

Centralizes magic numbers and default parameters used throughout the codebase
to improve maintainability and scientific reproducibility.
"""

# ============================================================================
# Random Seeds
# ============================================================================

DEFAULT_SEED = 43
"""Default random seed for reproducibility."""

SEED_BATCH_EFFECT = 43
"""Default seed for batch effect generation."""

# ============================================================================
# Fuzzy C-Means Parameters
# ============================================================================

FUZZY_CMEANS_M = 2.0
"""
Fuzziness parameter for fuzzy c-means clustering.

Controls the degree of fuzzy overlap between clusters:
- m = 1: Hard clustering (not allowed, causes division by zero)
- m = 2: Moderate fuzziness (standard choice, balanced)
- m → ∞: Maximum fuzziness (all memberships approach 1/k)

Typical range: 1.5 to 3.0
Default: 2.0 (most common in literature)
"""

FUZZY_CMEANS_EPS = 1e-10
"""
Epsilon value for numerical stability in fuzzy c-means.

Added to distances to prevent division by zero when:
- Cells are exactly at centroid locations
- All distances are very small
"""

# ============================================================================
# Sigmoid Transformation Parameters
# ============================================================================

SIGMOID_STEEPNESS = 0.5
"""
Steepness parameter (a) for sigmoid transformation of weights.

Formula: f(x) = 1 / (1 + exp(-a*x))

Controls how sharply weights transition from low to high:
- a = 0.1: Gentle sigmoid (gradual transition)
- a = 0.5: Moderate sigmoid (standard choice)
- a = 1.0: Steep sigmoid (sharp transition)
- a = 5.0: Very steep sigmoid (almost step function)

Higher values create more distinct cluster boundaries.
Default: 0.5
"""

# ============================================================================
# Probability Scaling Parameters
# ============================================================================

NEUTRAL_PROBABILITY = 0.5
"""
Neutral baseline probability for non-enriched cells.

In differential abundance analysis:
- prob > 0.5: Enriched in condition
- prob = 0.5: Neutral (no enrichment)
- prob < 0.5: Depleted (not used in current pipeline)
"""

MIN_ENRICHMENT = 0.5
"""Minimum allowed enrichment probability (neutral)."""

MAX_ENRICHMENT = 0.95
"""
Maximum typical enrichment probability.

While 1.0 is theoretically possible, 0.95 is used as a practical maximum
to avoid extreme certainty in label assignment.
"""

# ============================================================================
# Batch Effect Parameters
# ============================================================================

DEFAULT_BATCH_SD = 0.5
"""
Default batch effect strength as fraction of data variance.

Batch effects are added as:
batch_effect ~ Normal(0, batch_sd * sqrt(sum of component variances))

- batch_sd = 0: No batch effects
- batch_sd = 0.5: Moderate batch effects (default)
- batch_sd = 1.0: Strong batch effects (equal to data std)
- batch_sd = 2.0: Very strong batch effects

Typical range for benchmarking: 0.0 to 4.0
"""

# ============================================================================
# Pipeline Defaults
# ============================================================================

DEFAULT_N_CONDITIONS = 2
"""Number of experimental conditions (typically control vs treatment)."""

DEFAULT_N_REPLICATES = 3
"""Number of replicates per condition."""

DEFAULT_N_BATCHES = 2
"""Number of technical batches."""

# ============================================================================
# True Label Assignment Parameters
# ============================================================================

QUANTILE_LABEL_THRESHOLD_PCT = 10
"""
Percentage offset from enrichment level for true label assignment in cluster datasets.

Used to create da_lower and da_upper thresholds:
- If pop_enr > 0.5: da_upper = pop_enr - (pop_enr * threshold / 100)
- If pop_enr < 0.5: da_lower = pop_enr + (pop_enr * threshold / 100)

Default: 10 (i.e., 10% offset)
"""

# ============================================================================
# Numerical Tolerances
# ============================================================================

SCALE_ZERO_TOLERANCE = 1e-10
"""
Tolerance for determining if array has zero range in scale() function.

If max(x) - min(x) < tolerance, array is considered constant.
"""

PROBABILITY_SUM_TOLERANCE = 1e-6
"""
Tolerance for probability normalization checks.

Probabilities should sum to 1.0 within this tolerance.
"""

# ============================================================================
# Validation Ranges
# ============================================================================

VALID_ENRICHMENT_RANGE = (0.0, 1.0)
"""Valid range for enrichment probabilities (exclusive of endpoints)."""

VALID_M_RANGE = (1.0, float('inf'))
"""Valid range for fuzzy c-means m parameter (m must be > 1)."""

VALID_BATCH_SD_RANGE = (0.0, 10.0)
"""
Valid range for batch effect standard deviation.

While mathematically any positive value is valid, values > 10 are unrealistic
for benchmarking purposes (would completely swamp biological signal).
"""

# ============================================================================
# File Format Parameters
# ============================================================================

FLOAT_PRECISION = 16
"""Decimal places to save in CSV files for floating point values."""

# ============================================================================
# Helper Functions
# ============================================================================

def validate_enrichment(enr):
    """
    Validate enrichment probability value.

    Parameters:
    -----------
    enr : float
        Enrichment probability to validate

    Raises:
    -------
    ValueError
        If enrichment is outside valid range

    Returns:
    --------
    bool
        True if valid
    """
    if not (VALID_ENRICHMENT_RANGE[0] < enr < VALID_ENRICHMENT_RANGE[1]):
        raise ValueError(
            f"Enrichment must be in range {VALID_ENRICHMENT_RANGE}, got {enr}"
        )
    return True


def validate_m_parameter(m):
    """
    Validate fuzzy c-means m parameter.

    Parameters:
    -----------
    m : float
        Fuzziness parameter to validate

    Raises:
    -------
    ValueError
        If m <= 1 (invalid for fuzzy c-means)

    Returns:
    --------
    bool
        True if valid
    """
    if m <= VALID_M_RANGE[0]:
        raise ValueError(
            f"Fuzzy c-means parameter m must be > {VALID_M_RANGE[0]}, got {m}"
        )
    return True


def validate_batch_sd(batch_sd):
    """
    Validate batch effect standard deviation.

    Parameters:
    -----------
    batch_sd : float
        Batch effect strength to validate

    Raises:
    -------
    ValueError
        If batch_sd is negative or unrealistically large

    Returns:
    --------
    bool
        True if valid
    """
    if not (VALID_BATCH_SD_RANGE[0] <= batch_sd <= VALID_BATCH_SD_RANGE[1]):
        raise ValueError(
            f"Batch SD should be in range {VALID_BATCH_SD_RANGE}, got {batch_sd}"
        )
    return True
