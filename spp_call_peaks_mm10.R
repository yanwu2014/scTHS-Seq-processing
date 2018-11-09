#####################
# @desc No alts
# human visual cortex 
# (data-raw5)
# 
# @date Feburary 28, 2017
# @author Jean Fan
######################
# Modified by EDuong 1/31/2018
######################
# First call speaks using SPP
######################
#view all commands in spp: ls("package:spp")

require(spp)
require(parallel)
require(methods)

#### Load helper functions

# http://jef.works/blog/2016/12/06/mapping-snps-and-peaks-to-genes-in-R/
#' Convert from ranges to GRanges
#'
#' @param df Dataframe with columns as sequence name, start, and end
#'
#' @returns GRanges version
#'
range2GRanges <- function(df) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = df[,2], end = df[,3])
  )
  return(gr)
}


pos2GRanges <- function(chr, pos) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
  return(gr)
}


## Establish valid chromosmes
chrl <- paste0("chr", c(1:19, "X"))
names(chrl) <- chrl

#### Pooling data

## Get input parameters
options <- commandArgs(trailingOnly = T)
output.prefix <- options[[1]]
n.cores <- as.integer(options[[2]])
min.reads <- as.integer(options[[3]])
input.dirs <- options[4:length(options)]

## Format bam files
bamfiles <- unlist(lapply(input.dirs, function(input.path) {
  input.path <- gsub("/", "", input.path)
  print(input.path)
  
  bamfiles <- list.files(input.path, pattern = ".unique*bam$")
  names(bamfiles) <- gsub(".bam", "", bamfiles)
  print(paste0("Total number of files: ", length(bamfiles)))
  
  ## Count unique reads per cell
  reads <- unlist(mclapply(bamfiles, function(fname) {
    tryCatch({
      x <- read.bam.tags(paste(input.path, fname, sep = "/"))$tags
      sum(unlist(lapply(x[chrl], length)))
    }, error = function(r) return(NA)) ## empty files
  }, mc.cores = n.cores))
  names(reads) <- bamfiles
  reads <- na.omit(reads)
  
  ## Plot distribution of unique reads per cell
  pdf(paste0(input.path, "_unique_read_dist.pdf"), width = 5, height = 5)
  hist(log10(reads + 1), breaks = 50)
  abline(v = log10(min.reads), col = "red")
  dev.off()
  
  ## Only use files with more than the minimum number of reads
  bamfiles.use <- paste(input.path, names(which(log10(reads + 1) > log10(min.reads))), sep = "/")
  print(paste0("Total number of files: ", length(bamfiles.use)))
  
  bamfiles.use
}))
write(bamfiles, file = paste0(output.prefix, ".bamfiles.pass_filter.txt"), sep = "\n")

## Pool reads
bamdata <- mclapply(bamfiles, function(fname) {
  read.bam.tags(fname)$tags
}, mc.cores = n.cores)

## Pooled data. $tags will represent Tn5 cut positions on each chromosome (sign carries the strand information)
## $cell will give the id (integer) of a cell giving the read; both will be sorted in the chromosome coordinate order
pdata <- mclapply(chrl, function(chr) {
  list(tags = unlist(lapply(1:length(bamdata), function(i) na.omit(bamdata[[i]][[chr]]))), 
       cells = unlist(lapply(1:length(bamdata), function(i) rep(i, length(na.omit(bamdata[[i]][[chr]]))))))
}, mc.cores = n.cores)

## Sort by cut position
pdata <- mclapply(pdata, function(d) { 
  co <- order(abs(d$tags), decreasing = F)
  return(list(tags = d$tags[co], cells = d$cells[co]))
}, mc.cores = n.cores)

#### Remove reads in masked regions if necessary
library(BSgenome.Mmusculus.UCSC.mm10.masked)
genome <- BSgenome.Mmusculus.UCSC.mm10.masked

rl <- 100
pdata.filtered <- lapply(chrl, function(chr) {
  pdatachr <- pdata[[chr]] ## unfiltered
  ## expand mask region to account for read length
  mask <- GenomicRanges::GRanges(seqnames = chr, 
                                 IRanges(start = masks(genome[[chr]])$RM@start - rl,
                                         end = masks(genome[[chr]])$RM@start + masks(genome[[chr]])$RM@width + rl))
  p <- pos2GRanges(chr, abs(pdatachr$tags))
  overlap <- GenomicRanges::findOverlaps(p, mask, ignore.strand=T)

  ## remove these bad reads
  sd <- setdiff(1:length(pdatachr$tags), slot(overlap, 'from'))
  pdatachr$tags <- pdatachr$tags[sd]
  pdatachr$cells <- pdatachr$cells[sd]

  return(pdatachr)
})
names(pdata.filtered) <- chrl

##### Call Peaks

## Peak calling parameters
bandwidth <- 500
step <- 100
thr <- 5
span <- 10
fdr <- 1e-8


## Call peaks
smoothed.density <- lapply(lapply(pdata.filtered, function(d) abs(d$tags)), function(d) {
  tc <- window.tag.count(d, window.size = bandwidth, window.step = step)
  x <- seq(tc$x[1], tc$x[2], by = tc$step)
  y <- tc$y
  data.frame("x" = x, "y" = y)
})    
names(smoothed.density) <- chrl

## this will calculate index positions of all local maxima on the density vector $y
peaks.all <- lapply(chrl, function(chr) {
  peak.indices <- peaks.c(smoothed.density[[chr]]$y, thr = thr, max.span = span)
  df <- data.frame(peak.position = smoothed.density[[chr]]$x[peak.indices], 
                   peak.magnitude = smoothed.density[[chr]]$y[peak.indices])    
  return(df)
})
names(peaks.all) <- chrl
print(paste0("Number of raw peaks: ", sum(unlist(lapply(peaks.all, nrow)))))

## Filter peaks
peaks.filtered <- mclapply(chrl, function(chr) {
  test.x <- smoothed.density[[chr]]$x
  test.y <- smoothed.density[[chr]]$y
  
  ## init peaks
  df <- peaks.all[[chr]]
  fn <- ecdf(df$peak.magnitude)
  
  ## randomize
  set.seed(0)
  tags <- pdata.filtered[[chr]]$tags
  shuffle <- sort(runif(length(tags), min = range(tags)[1], max = range(tags)[2]), decreasing = TRUE)
  
  shuffle.smooth <- lapply(list(chr = abs(shuffle)), function(d) {
    tc <- window.tag.count(d, window.size = bandwidth, window.step = step)
    x <- seq(tc$x[1], tc$x[2], by = tc$step)
    y <- tc$y
    data.frame('x'= x,'y'= y)
  })
  peak.indices.shuffle <- peaks.c(shuffle.smooth[[1]]$y, thr = thr, max.span = span);
  df.shuffle <- data.frame(peak.position = shuffle.smooth[[1]]$x[peak.indices.shuffle], peak.magnitude = shuffle.smooth[[1]]$y[peak.indices.shuffle] );
  fn.shuffle <- ecdf(df.shuffle$peak.magnitude)
  
  t <- knots(fn.shuffle)
  
  V <- (1 - fn.shuffle(t))*nrow(df.shuffle)
  S <- (1 - fn(t))*nrow(df) 
  
  pseudo <- 1e-6
  fdrs <- (V + pseudo) / (V + S + pseudo)
  
  ## optimal threshold
  above <- fdrs > fdr
  intersect.points <- which(diff(above) != 0)
  threshold <- t[intersect.points][1]
  
  table(df$peak.magnitude > threshold)
  table(df.shuffle$peak.magnitude > threshold)
  
  df <- peaks.all[[chr]]
  df.t <- df[df$peak.magnitude > threshold,]
  return(df.t)
}, mc.cores = length(chrl))
names(peaks.filtered) <- chrl
print(paste0("Number of filtered peaks: ", sum(unlist(lapply(peaks.filtered, function(x) nrow(na.omit(x)))))))


## spp results
peaks.spp <- do.call(rbind, lapply(chrl, function(chr) {
  df <- peaks.filtered[[chr]]
  data.frame(chr, df[,1] - bandwidth/2, df[,1] + bandwidth/2)
}))
peaks.spp <- peaks.spp[complete.cases(peaks.spp),]
peaks.spp <- range2GRanges(peaks.spp)
names(peaks.spp) <- paste(peaks.spp)

## write out bed file
gr <- peaks.spp
bed.df <- data.frame(seqnames = seqnames(gr),
                     starts = start(gr) - 1,
                     ends = end(gr),
                     names = c(rep(".", length(gr))),
                     scores = c(rep(".", length(gr))),
                     strands = strand(gr))

write.table(bed.df, file = paste0(output.prefix, "_peaks_spp", ".bed"), 
            quote = F, sep = "\t", row.names = F, col.names = F)

## write out gtf file
peak.names <- paste(bed.df$seqnames, bed.df$starts, bed.df$ends, sep = "_")
gtf.df <- data.frame(chr = seqnames(gr), 
                     source = rep("spp", length(gr)),
                     feature = rep("peak", length(gr)),
                     starts = start(gr) - 1,
                     ends = end(gr),
                     scores = rep(".", length(gr)),
                     strands = rep("+", length(gr)),
                     frames = rep(".", length(gr)),
                     attr = paste0('peak_id \"', peak.names, '\";'))

write.table(gtf.df, file = paste0(output.prefix, "_peaks_spp", ".gtf"), 
            quote = F, sep = "\t", row.names = F, col.names = F)