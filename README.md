# eTFCE

**Exact Threshold-Free Cluster Enhancement via Fast Cluster Retrieval (eTFCE)** is a statistical inference procedure that integrates exact TFCE computation, rather than the approximate methods used in common implementations, with an efficient cluster retrieval algorithm. This enables fast nonparametric inference as well as simultaneous cluster-level inference.

## Installation

1. Install R: https://cran.r-project.org/
   
   Install RStudio (optional): https://posit.co/download/rstudio-desktop/

2. Install the R package `devtools`:
    ```r
    install.packages("devtools")
    ```

3. Install eTFCE from GitHub:
    ```r
    devtools::install_github("xuchen312/eTFCE")
    ```

## Examples

```r
library(eTFCE)

# Example 1: run eTFCE using automatically generated sign-flipping permutations
res <- eTFCE_perms(conts, mask)

# Example 2: run eTFCE using a user-provided permutation matrix (n_perm x n_subj)
sign_matrix <- read.table("TFCE_perm_tstat1.txt", header = FALSE)
sign_matrix <- unlist(t(sign_matrix), use.names = FALSE)
res <- eTFCE_perms(conts, mask, sign_matrix = sign_matrix, auto_signflip = FALSE)
```

## Input / Output

**Input**: 2 required arguments and 8 optional parameters

- `conts`: Path to 4D NIfTI image (nDim[1]*nDim[2]*nDim[3]*n_subj)
- `mask`: Path to 3D mask
- `n_perm`: Number of permutations to perform (default: 5000)
- `E`: TFCE extent parameter (default: 0.5)
- `H`: TFCE height parameter (default: 2.0)
- `t0`: Initial cluster-forming threshold (default: 3.1)
- `conn`: Connectivity criterion (e.g., 6, 18, 26; default: 26)
- `seed`: Random seed for reproducibility (default: 5000)
- `sign_matrix`: Vectorized permutation sign matrix (default: NULL). Used only when `auto_signflip = FALSE`.
- `auto_signflip`: Logical (default: TRUE). If `TRUE`, permutations are generated automatically via internal sign flipping. If `FALSE`, a user-provided `sign_matrix` must be supplied.

**Output**: a list of permutation-based null distributions and observed test statistics:

- `n_perm`: Number of permutations used.
- `maxTFCE`: Maximum TFCE statistic across permutations.
- `maxCZ`: Maximum cluster extent (CZ) across permutations.
- `maxCM`: Maximum cluster mass (CM) across permutations.
- `mTFCE`: Maximum TFCE value from the original (observed) data.
- `orgTFCE`: TFCE values from the original data.
- `orgCZ`: Cluster extent values from the original data.
- `orgCM`: Cluster mass values from the original data.

## Notes

- This repository currently provides an example for a one-sample t-test using sign flipping.
- The underlying computational framework is general and can be extended to other designs (e.g., two-sample tests or regression models). Extensions to other models can be provided upon request.

## References

Chen X., Weeda W., Nichols T.E. and Goeman J.J. (2026). eTFCE: Exact Threshold-Free Cluster Enhancement via Fast Cluster Retrieval. arXiv:2603.03004. doi:10.48550/arXiv.2603.03004.

Chen X., Goeman J.J., Krebs T.J., Meijer R.J. and Weeda W.D. (2023). Adaptive Cluster Thresholding with Spatial Activation Guarantees Using All-resolutions Inference. arXiv:2206.13587. doi:10.48550/arXiv.2206.13587.

Smith S.M. and Nichols T.E. (2009). Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. *NeuroImage*, 44(1), 83–98. doi:10.1016/j.neuroimage.2008.03.061

## Contact

If you are interested in applying this method to your own data or would like help adapting the code, please feel free to open an issue or get in touch. I would be happy to collaborate. :blush:

Email: [xuchen312@gmail.com](mailto:xuchen312@gmail.com)
