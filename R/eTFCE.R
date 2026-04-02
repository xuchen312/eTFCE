#' Run eTFCE permutation analysis using one-sample t-test
#'
#' This function performs nonparametric eTFCE inference using either
#' auto-generated sign-flipping permutations or a user-provided sign matrix.
#'
#' @name run_eTFCE
#'
#' @param conts Path to 4D NIfTI image
#' @param mask Path to 3D mask
#'
#' @param n_perm Number of permutations to perform (default: 5000)
#' @param E TFCE extent parameter (default: 0.5)
#' @param H TFCE height parameter (default: 2.0)
#' @param t0 Initial cluster-forming threshold (default: 3.1)
#' @param conn Connectivity criterion (e.g., 6, 18, 26; default: 26)
#' @param seed Random seed for reproducibility (default: 5000)
#' @param sign_matrix Optional permutation sign matrix (vectorized).
#'        (default: NULL)
#'        Used only when `auto_signflip = FALSE`.
#'        If provided while `auto_signflip = FALSE`, this matrix will be
#'        used for permutations instead of generating sign-flipping
#'        permutations.
#'        The matrix must have dimensions (n_perm x n_subjects)
#'        in column-major vectorized order.
#' @param auto_signflip Logical. (default: TRUE)
#'        If TRUE, permutations are generated automatically via sign flipping.
#'        If FALSE, a user-provided `sign_matrix` must be supplied.
#'
#' @return eTFCE inference results returned from the C++ backend, including
#'         permutation-based null distributions and observed test statistics.
#' \describe{
#'   \item{n_perm}{Number of permutations used.}
#'   \item{maxTFCE}{Maximum TFCE statistic across permutations.}
#'   \item{maxCZ}{Maximum cluster extent (CZ) across permutations.}
#'   \item{maxCM}{Maximum cluster mass (CM) across permutations.}
#'   \item{mTFCE}{Maximum TFCE value from the original (observed) data.}
#'   \item{orgTFCE}{TFCE values from the original data.}
#'   \item{orgCZ}{Cluster extent values from the original data.}
#'   \item{orgCM}{Cluster mass values from the original data.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # =========================================================
#' # Example 1: Auto sign-flipping permutations
#' # =========================================================
#'
#' res <- eTFCE_perms(conts, mask, n_perm = 1000, seed = 1000)
#'
#' # =========================================================
#' # Example 2: Use external sign matrix (e.g., FSL randomise)
#' # =========================================================
#'
#' sign_matrix <- read.table("TFCE_perm_tstat1.txt", header = FALSE)
#' sign_matrix <- unlist(sign_matrix, use.names = FALSE)
#'
#' res <- eTFCE_perms("subs_4D.nii.gz", "mask.nii.gz", n_perm = 1000,
#'                    sign_matrix = sign_matrix, auto_signflip = FALSE)
#'
#' }
#' 
#' 

run_eTFCE <- function(conts, mask, n_perm=5000, E=0.5, H=2.0, t0=3.1, conn=26, 
                      seed=5000, sign_matrix=NULL, auto_signflip=TRUE) {
  
  # load data & mask
  conts <- RNifti::readNifti(conts)
  mask <- RNifti::readNifti(mask)
  
  # extract dimensions & number of subjects
  dims4D <- dim(conts)
  dims <- dims4D[1:3]
  n <- dims4D[4]
  
  # reshape data for C++ function input
  conts <- c(conts)
  mask <- as.integer(c(mask))
  
  # run eTFCE inference
  eTFCE_perms(conts, mask, dims, n, n_perm=n_perm, E=E, H=H, t0=t0, conn=conn, 
              seed=seed, sign_matrix=sign_matrix, auto_signflip=auto_signflip)
}
