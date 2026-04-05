#' Run eTFCE permutation analysis using one-sample t-test
#'
#' This function performs nonparametric eTFCE inference using either
#' auto-generated sign-flipping permutations or a user-provided sign matrix.
#'
#' @name run_eTFCE
#'
#' @param conts Path to a 4D NIfTI image
#' @param mask Path to a 3D mask
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
#' @return A list of eTFCE inference results returned from the C++ backend, 
#'         including permutation-based null distributions and observed test 
#'         statistics.
#' \describe{
#'   \item{n_perm}{Number of permutations used (scalar).}
#'   \item{maxTFCE}{Maximum TFCE statistic across permutations (vector).}
#'   \item{maxCZ}{Maximum cluster extent (CZ) across permutations (vector).}
#'   \item{maxCM}{Maximum cluster mass (CM) across permutations (vector).}
#'   \item{mTFCE}{Maximum TFCE value from the original (observed) data (scalar).}
#'   \item{orgCZ}{Cluster extent values from the original data (vector).}
#'   \item{orgCM}{Cluster mass values from the original data (vector).}
#'   \item{orgTFCE}{TFCE values from the original data (vector/image).}
#'   \item{orgT}{t-statistic values from the original data (vector/image).}
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
#' res <- run_eTFCE(conts, mask, n_perm = 1000, seed = 1000)
#'
#' # =========================================================
#' # Example 2: Use external sign matrix (e.g., FSL randomise)
#' # =========================================================
#'
#' # sign matrix (n_perm × n_subj), saved as a text file 
#' # (output from FSL randomise).
#' sign_matrix <- read.table("TFCE_perm_tstat1.txt", header = FALSE)
#' sign_matrix <- unlist(t(sign_matrix), use.names = FALSE)
#'
#' res <- run_eTFCE(conts, mask, n_perm = 1000,
#'                  sign_matrix = sign_matrix, auto_signflip = FALSE)
#'
#' }
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

# extract null dist
get_null_dist <- function(x, maxx) {
  if (abs(x[1] - maxx) < 1e-12) return(x[-1])
  return(x)
}

#' Write eTFCE results to NIfTI files and generate cluster-level result table
#'
#' This function writes voxel-wise statistical maps (t-statistic, TFCE,
#' and TFCE-based corrected p-values) to NIfTI files, and computes
#' cluster-level inference results.
#' 
#' @name write_eTFCE
#'
#' @param res A list of eTFCE inference results returned by \code{run_eTFCE}.
#' @param mask Path to a 3D mask.
#' @param file Output file prefix (character). Suffixes will be appended
#'        automatically.
#'
#' @return A data.frame summarizing cluster-level inference results.
#' \describe{
#'   \item{Cluster}{Cluster label (e.g., cl1, cl2, ...).}
#'   \item{CZ}{Cluster extent from the original data.}
#'   \item{p_CZ}{FWER-corrected p-values for cluster extent.}
#'   \item{CM}{Cluster mass from the original data.}
#'   \item{p_CM}{FWER-corrected p-values for cluster mass.}
#' }
#'
#' @details
#' The function outputs the following NIfTI files:
#' \itemize{
#'   \item \code{*_tstat1.nii.gz}: voxel-wise t-statistic map
#'   \item \code{*_tfce.nii.gz}: voxel-wise TFCE map
#'   \item \code{*_tfce_corrp.nii.gz}: TFCE-based FWER-corrected p-value map
#' }
#'
#' Voxel-wise p-values are computed using permutation-based null distributions
#' of the maximum TFCE statistic.
#'
#' @export
#' 
write_eTFCE <- function(res, mask, file) {
  
  mask <- RNifti::readNifti(mask)
  dims <- dim(mask)
  
  # reshape & output statistic maps
  orgT    <- array(res$orgT,    dim=dims)
  orgTFCE <- array(res$orgTFCE, dim=dims)
  
  RNifti::writeNifti(orgT,    file=paste0(file,"_tstat1.nii.gz"), template=mask)
  RNifti::writeNifti(orgTFCE, file=paste0(file,"_tfce.nii.gz"),   template=mask)
  
  # ===============================
  # TFCE-based FWER-corrected p-map
  # ===============================
  
  # vectorized mask & TFCE map
  mask_vec <- as.logical(c(mask))
  tfce_vec <- res$orgTFCE
  
  # check if original result is within max null distribution
  null_TFCE <- get_null_dist(res$maxTFCE, res$mTFCE)
  n_null <- length(null_TFCE)
  
  # compute FWER-corrected p-values
  pvalTFCE <- rep(1, length(tfce_vec))
  pvalTFCE[mask_vec] <- sapply(
    tfce_vec[mask_vec],
    function(x) (sum(null_TFCE >= x) + 1) / (n_null + 1)
  )
  pvalTFCE <- array(pvalTFCE, dim=dims)
  
  RNifti::writeNifti(pvalTFCE, paste0(file, "_tfce_corrp.nii.gz"), template=mask)
  
  # =======================
  # cluster-level inference
  # =======================
  nclus <- length(res$orgCZ)
  
  null_CZ <- get_null_dist(res$maxCZ, max(res$orgCZ))
  null_CM <- get_null_dist(res$maxCM, max(res$orgCM))
  n_null <- length(null_CZ)
  
  pvalCZ <- sapply(res$orgCZ, function(x) (sum(null_CZ >= x)+1) / (n_null+1))
  pvalCM <- sapply(res$orgCM, function(x) (sum(null_CM >= x)+1) / (n_null+1))
  
  # ===================
  # output result table
  # ===================
  df <- data.frame(Cluster=paste0("cl", seq_len(nclus)),
                   CZ=res$orgCZ, p_CZ=pvalCZ,
                   CM=res$orgCM, p_CM=pvalCM)
  df <- df[order(-df$CZ),]
  row.names(df) <- NULL
  
  return(df)
}
