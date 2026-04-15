# eTFCE

**Exact Threshold-Free Cluster Enhancement via Fast Cluster Retrieval (eTFCE)** is a statistical inference approach for neuroimaging data that combines exact TFCE computation with an efficient cluster retrieval algorithm.

Unlike standard TFCE implementations that rely on approximations, eTFCE performs exact TFCE-based permutation inference, enabling both voxel-wise and cluster-level nonparametric inference in a unified computational framework.

The toolbox also provides utilities for statistical map reconstruction and inference output generation, including TFCE-corrected p-value maps and cluster-level summary tables.

## Installation

1. Install R: https://cran.r-project.org/
   
   Install RStudio (optional): https://posit.co/download/rstudio-desktop/

   Install the R package `devtools`:
    ```r
    install.packages("devtools")
    ```

2. Install eTFCE from GitHub:
    ```r
    devtools::install_github("xuchen312/eTFCE")
    ```

## Usage

### Example 1: Automatic sign-flipping permutations

```r
library(eTFCE)

res <- eTFCE_perms(conts, mask, n_perm = 1000)
```

### Example 2: User-provided permutation matrix (nSubj x nPerm)

```r
# FSL randomise outputs nPerm x nSubj permutation matrix
sign_matrix <- read.table("TFCE_perm_tstat1.txt", header = FALSE)
# transpose and vectorize
sign_matrix <- unlist(t(sign_matrix), use.names = FALSE)

res <- eTFCE_perms(
conts, mask,
sign_matrix = sign_matrix,
auto_signflip = FALSE)
```

### Example 3: Write results

After running inference, you can generate statistical maps and cluster tables:

```r
df <- write_eTFCE(res, mask, file = "results/eTFCE_output")
```

This produces:

- `results/eTFCE_output_tstat.nii.gz` - Voxel-wise t-statistic map
- `results/eTFCE_output_tfce.nii.gz` - Voxel-wise TFCE map
- `results/eTFCE_output_tfce_corrp.nii.gz` - TFCE-corrected p-value map
- Cluster-level inference result table (printed to console; also returned as data frame)

To save the cluster table as a CSV file:

```r
write.csv(df, "clus.csv", row.names = FALSE)
```

`write_eTFCE()` is designed to be robust to inconsistent permutation outputs (i.e., whether observed statistics are included in permutation-based output vectors or not).

## Input Parameters

### Required inputs

- `conts`: Path to 4D NIfTI image (nDim[1] x nDim[2] x nDim[3] x n_subj), which can be generated from multiple 3D NIfTI images using FSL, e.g.,

   ```r
   fslmerge -t output_4D.nii.gz sub*.nii.gz
   ```

- `mask`: Path to 3D mask

### Optional parameters

- `n_perm`: Number of permutations (default: 5000)
  
- `E`: TFCE extent parameter (default: 0.5)
  
- `H`: TFCE height parameter (default: 2.0)
  
- `t0`: Initial cluster-forming threshold (default: 3.1)
  
- `conn`: Voxel connectivity (6 / 18 / 26; default: 26)
  
- `seed`: Random seed for reproducibility (default: 5000)
  
- `sign_matrix`: Vectorized User-supplied (nSubj x nPerm) permutation sign matrix (default: NULL)
  
- `auto_signflip`: If `TRUE`, permutations are generated internally via sign-flipping; if `FALSE`, a user-provided `sign_matrix` must be supplied (default: TRUE)

## Output of `eTFCE_perms`

A list of permutation-based null distributions and observed statistics:

### Permutation-based statistics

- `n_perm`: Number of permutations used

- `maxTFCE`: Maximum TFCE statistic across permutations
  
- `maxCZ`: Maximum cluster extent (CZ) across permutations
  
- `maxCM`: Maximum cluster mass (CM) across permutations

### Observed statistics
  
- `mTFCE`: Maximum TFCE value from the original (observed) data
  
- `orgCZ`: Cluster extent values (observed)
  
- `orgCM`: Cluster mass values (observed)
  
- `orgTFCE`: Voxel-wise TFCE map (observed)
  
- `orgT`: Voxel-wise t-statistic map (observed)

The returned list includes both permutation-based null distributions and observed statistics. The exact structure is handled internally and is compatible with `write_eTFCE()` for post-processing.

## Post-processing (`write_eTFCE`)

The function `write_eTFCE()` provides:

1. Voxel-wise outputs

   - t-statistic map
   - TFCE map
   - TFCE-corrected p-value map (FWER)
  
2. Cluster-level inference

   - cluster extent (CZ)
   - cluster mass (CM)
   - permutation-based p-values
  
3. Output table

   Returned table:
   
   Cluster | CZ | p_CZ | CM | p_CM 
   --------|----|------|----|------
   ...|...|...|...|...

   Clusters are sorted by decreasing cluster extent (CZ).

## Notes

- The current implementation supports **one-sample t-test via sign-flipping permutations**
- The framework is general and can be extended to:
      - two-sample tests
      - regression models
      - GLM-based designs
      - etc.

## References

Chen X., Weeda W., Nichols T.E. and Goeman J.J. (2026). eTFCE: Exact Threshold-Free Cluster Enhancement via Fast Cluster Retrieval. arXiv:2603.03004. doi:10.48550/arXiv.2603.03004.

Chen X., Goeman J.J., Krebs T.J., Meijer R.J. and Weeda W.D. (2023). Adaptive Cluster Thresholding with Spatial Activation Guarantees Using All-resolutions Inference. arXiv:2206.13587. doi:10.48550/arXiv.2206.13587.

Smith S.M. and Nichols T.E. (2009). Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. *NeuroImage*, 44(1), 83–98. doi:10.1016/j.neuroimage.2008.03.061

## Contact

If you are interested in applying this method to your own data or extending it to other study designs, please open an issue or contact the author. Contributions and collaborations are welcome. :blush:

Email: [xuchen312@gmail.com](mailto:xuchen312@gmail.com)
