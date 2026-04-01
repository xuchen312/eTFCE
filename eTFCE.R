library(RNifti)
source("C:/Users/xc23776/Downloads/eTFCE/eTFCE_funs.R")

# ---------- [0] Preparation (if needed)
filepath <- "C:/Users/xc23776/Downloads/fMRIdata/data-raw/AuditoryData"
mask <- readNifti("C:/Users/xc23776/Downloads/fMRIdata/data-raw/mask.nii.gz")
mask <- mask!=0
nii_files <- list.files(
  path=file.path(filepath),
  pattern="\\.nii(\\.gz)?$",
  full.names=TRUE,
  recursive=TRUE
)
n <- length(nii_files)
dims <- RNifti::niftiHeader(nii_files[1])$dim[2:4]
idt <- which(mask)

# read in all images
imgs <- lapply(nii_files, RNifti::readNifti)
copes <- array(unlist(imgs, use.names=FALSE), dim=c(dims,n))

# compute original t-map
avg_tmap <- apply(copes, c(1,2,3), mean)
std_tmap <- apply(copes, c(1,2,3), sd)
tmap <- sqrt(n) * avg_tmap / std_tmap
tmap[std_tmap==0] <- 0
mask <- tmap!=0
RNifti::writeNifti(tmap, "tmap.nii.gz", template=copes[,,,1])

# ---------- [1] Compute original TFCE scores & cluster statistics
etfce_org <- eTFCEorg(tmap, mask, conn=26)
etfceTFCE_org <- etfce_org$orgTFCE
etfceCE_org <- etfce_org$orgCE
etfceCM_org <- etfce_org$orgCM

# -------- [2] Compute max(T) null distributions
n_perm <- 1000
etfce_perm <- eTFCE_perms(conts=copes, mask=mask, n_perm=n_perm)

# -------- [3] Compute FWER-corrected p-values

# ------ For TFCE
tfce <- etfceTFCE_org[idt]
pval_tfce <- (rowSums(outer(tfce, etfce_perm$maxTFCE, "<=")) + 1) / (etfce_perm$n_perm + 1)
# map p-values back to the 3D space
pval_eTFCE <- array(1, dims)
pval_eTFCE[idt] <- pval_tfce
RNifti::writeNifti(pval_eTFCE, "pval_etfce.nii.gz", template=tmap)

# ------ For cluster extent & cluster mass
pval_ce <- sapply(etfceCE_org, function(xi) sum(xi<=etfce_perm$maxCZ)+1) / (etfce_perm$n_perm+1)
pval_cm <- sapply(etfceCM_org, function(xi) sum(xi<=etfce_perm$maxCM)+1) / (etfce_perm$n_perm+1)




