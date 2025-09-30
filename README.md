# LocalizedSampling
Localized Sampling Algorithm

### Overview
Implements a localized sampling workflow that:
- partitions the numeric feature space,
- estimates within‑partition class heterogeneity (Shannon entropy), and
- samples a final set that prioritizes heterogeneous partitions while still covering homogeneous ones, using systematic samplers.

Dependencies: `dplyr`, `tidyr`, `ggplot2`, `tibble`, `purrr`.

### Functions

#### `localized_sampling(data, y_var, final_sample_size = NULL, final_sample_prop = NULL, n_partitions = 6, heterogeneous_prop = 0.8, return_idx = FALSE)`
- **Purpose**: End‑to‑end pipeline. Creates partitions, computes class proportions and entropy by partition, splits the target sample between heterogeneous and homogeneous partitions, and draws systematic samples from each. Ensures at least one instance of every class is present in the final sample.
- **Arguments**
  - **data**: data.frame/tibble with numeric features and a class column.
  - **y_var**: name of the class/label column in `data`.
  - **final_sample_size**: desired total number of rows in the final sample (exclusive with `final_sample_prop`).
  - **final_sample_prop**: proportion of the full data to sample (0 < p ≤ 1) if `final_sample_size` is not given.
  - **n_partitions**: number of bins per numeric feature used to form partitions.
  - **heterogeneous_prop**: fraction of the final sample allocated to partitions with entropy > 0.
  - **return_idx**: currently returns sampled rows; indices not returned.
- **Returns**: tibble with sampled rows from both heterogeneous and homogeneous partitions.

#### `create_partitions(df, n_partitions = 4)`
- **Purpose**: Cuts each numeric feature into `n_partitions` bins and encodes the multi‑dimensional bin as a single integer partition id.
- **Arguments**
  - **df**: input data.frame/tibble.
  - **n_partitions**: number of bins per numeric feature.
- **Returns**: factor vector of partition ids aligned with `nrow(df)`.

#### `class_props(df, part = partition, y = y, wide = FALSE)`
- **Purpose**: Computes per‑partition class counts and proportions, preserving zero‑count classes via `tidyr::complete`. Optionally pivots to wide format of class proportions.
- **Arguments**
  - **df**: input data including partition and class columns.
  - **part**: tidyselect for partition column (default `partition`).
  - **y**: tidyselect for class column (default `y`).
  - **wide**: if `TRUE`, returns a wide matrix of class proportions by partition.
- **Returns**: long tibble with counts and proportions, or a wide tibble if `wide = TRUE`.

#### `class_entropy(probs)`
- **Purpose**: Computes Shannon entropy of a probability vector, ignoring zeros.
- **Arguments**
  - **probs**: numeric vector of class probabilities that sum to 1.
- **Returns**: numeric entropy value in bits.

#### `sys_o2(data, y_var, n_total = NULL, prop = NULL, order_by = NULL, decreasing = FALSE, within = c("random", "given"), seed = NULL, return_idx = FALSE)`
- **Purpose**: Systematic sampling over within‑class order. Allocates a sample size to each class by proportional allocation with largest remainders, then selects entries systematically within each class.
- **Arguments**
  - **data**: data.frame/tibble.
  - **y_var**: class column name.
  - **n_total**: desired total sample size (exclusive with `prop`).
  - **prop**: global sampling proportion if `n_total` not provided.
  - **order_by**: optional column to order within each class before systematic selection.
  - **decreasing**: if `TRUE`, sorts `order_by` decreasing.
  - **within**: "random" to shuffle within class when `order_by` is `NULL`; "given" keeps input order.
  - **seed**: random seed for reproducibility when randomization occurs.
  - **return_idx**: if `TRUE`, returns row indices; otherwise returns sampled rows with metadata columns (e.g., `.rank_sys`).
- **Returns**: sampled rows or indices; includes allocation metadata when returning rows.

#### `sys_pc(data, y_var = NULL, n_total = NULL, prop = NULL, scale_vars = TRUE, seed = NULL, return_idx = FALSE)`
- **Purpose**: Systematic sampling from a single partition based on the angle around the partition centroid in a 1D/2D PCA score space, promoting spatial coverage.
- **Arguments**
  - **data**: data for one partition.
  - **y_var**: optional response/label column to exclude from features when building PCA.
  - **n_total**: desired number from this partition (exclusive with `prop`).
  - **prop**: proportion of rows to sample if `n_total` not given.
  - **scale_vars**: if `TRUE`, standardizes numeric features before PCA.
  - **seed**: random seed.
  - **return_idx**: if `TRUE`, returns row indices; else returns sampled rows with `.pc1`, `.pc2`, `.angle`, `.rank_sys`.
- **Returns**: sampled rows or indices from the given partition.

### Minimal usage
```r
# df: a data.frame with numeric features and a class column named "y"
set.seed(123)
sampled <- localized_sampling(
  data = df,
  y_var = "y",
  final_sample_prop = 0.2,
  n_partitions = 6,
  heterogeneous_prop = 0.8
)
```
