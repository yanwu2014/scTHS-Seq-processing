library(Rsubread)
library(Matrix)
library(methods)

## Command line arguments
options <- commandArgs(trailingOnly = T)
bamfiles.list <- options[[1]]
peaks.gtf.file <- options[[2]]
output.prefix <- options[[3]]
n.cores <- as.integer(options[[4]])

## Measure runtime
start.time <- proc.time()

## Count reads in peaks
max.files <- 5000
bamfiles.all <- scan(bamfiles.list, what = character(), sep = "\n")
bamfiles.list <- split(bamfiles.all, ceiling(seq_along(bamfiles.all)/max.files))

## Run read counting on each batch of 5000 cells
counts.list <- lapply(bamfiles.list, function(bamfiles) {
  ## Run feature counts
  res <- featureCounts(bamfiles, annot.ext = peaks.gtf.file, isGTFAnnotationFile = T,
                       GTF.featureType = "peak", GTF.attrType = "peak_id", ignoreDup = T,
                       nthreads = n.cores)
  
  ## Pull out metadata
  stats <- res$stat; rownames(stats) <- stats$Status; stats$Status <- NULL; stats <- as.matrix(stats);
  frac.in.peaks <- stats["Assigned",]/stats["Unassigned_NoFeatures",]
  
  ## Return counts matrix with metadata row
  as(rbind(frac.in.peaks, res$counts), "dgCMatrix")
})
counts <- do.call(cbind, counts.list)

## Print runtime
print("Time elapsed:")
print(proc.time() - start.time)
rm(counts.list); invisible(gc());

reads.in.peaks <- counts[1,]
counts <- counts[2:nrow(counts),]

save(counts, reads.in.peaks, file = paste0(output.prefix, ".counts.RData"))