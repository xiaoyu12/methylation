library(Gviz)
library(ggbio)

# Calculate methylation m values from methylBase object
calc_M_values <- function (meth) {
  t <- getData(meth)         # extract data part into a data frame
  nsamples <- (length(t) - 4) / 3
  meth.data <- t[, 1:4]
  for (i in 1:nsamples) {
    colname = paste0("m", i)
    # col 3i+3 for ith numC's, col 3i+4 for ith numT's 
    meth.data[, colname] <- log2((t[, 3*i+3] + 1)/(t[, 3*i+4] + 1)) # Add one to avoid dividing by 0
  }
  return(meth.data)
}
  
# Convert a region methylRaw list into a data.frame
# methylraw:  Large methylRawList 
# regions: data frame of genomic regions
methyl_to_data_frame <- function(methylraw, regions) {
  n <- length(methylraw)  # length of the list
  methyl_list <- as.list(1:n)
  for (i in 1:n) {
    t <- merge(getData(methylraw[[i]]), regions) # merge data frames according to "chr", "start", "end" and "strand"
    t <- t[!duplicated(t$ID), ]     # Remove duplicated ID's in the data frame
    t <- t[t$coverage >= 100, ]
    t$beta <- t$numCs / t$coverage
    t$m <- log2((t$numCs + 1) / (t$numTs + 1))  # Add one to avoid dividing by 0
    rownames(t) <- t$ID
    methyl_list[[i]] <- t
  }
  
  # Find shared IDs of data frames in methyl_list
  ids <- methyl_list[[1]]$ID
  for (i in 2:n) {
    ids <- intersect(ids, methyl_list[[i]]$ID)
  }
  # Get m values of the first sample
  data <- as.data.frame(methyl_list[[1]][as.character(ids), "m"])
  for (i in 2:n) {
    data <- cbind(data, methyl_list[[i]][as.character(ids), "m"])
  }
  colnames(data) <- paste0("m", 1:n)
  return(data)
}
  
# Get methylation values for gene annotated features: e.g. promoters, exons ...
getFeatureMethyl <- function(m.obj, g.features, s.names, lo.count=100) {
  mapping <- read.table("../v2/genbank_mapping.txt", row.names = 1)
  colnames(mapping) <- c("gene", "ID")
  nsamples <- length(m.obj)
  promoters <- as.data.frame(g.features)
  promoters$ID = mapping[promoters$name, "ID"]
  promoters = promoters[!duplicated(promoters$ID), ]    # remove duplicates due to alternative splicing
  colnames(promoters)[1] <- "chr"
  
  promobj <- regionCounts(m.obj, g.features)
  
  meth.prom <- list()
  for (i in 1:nsamples) {
    t <- merge(getData(promobj[[i]]), promoters, by.x=c("chr", "start", "end", "strand"), by.y=c("chr", "start", "end", "strand"))
    t <- t[!duplicated(t$ID), ]    # Remove duplicates due to alternative splicing
    t <- t[t$coverage >= lo.count, ]
    t$beta <- t$numCs / t$coverage
    t$m <- log2((t$numCs+1)/(t$numTs+1))
    rownames(t) <- t$ID
    meth.prom[[i]] <- t
  }
  
  ids <- meth.prom[[1]]$ID
  for (i in 2:nsamples) {
    ids <- intersect(ids, meth.prom[[i]]$ID)
  }
  mtx.coverage <- meth.prom[[1]][as.character(ids), "coverage"]
  mtx.m <- meth.prom[[1]][as.character(ids), "m"]
  mtx.beta <- meth.prom[[1]][as.character(ids), "beta"]
  for (i in 2:nsamples) {
    mtx.coverage <- cbind(mtx.coverage, meth.prom[[i]][as.character(ids), "coverage"])
    mtx.m <- cbind(mtx.m, meth.prom[[i]][as.character(ids), "m"])
    mtx.beta <- cbind(mtx.beta, meth.prom[[i]][as.character(ids), "beta"])
  }
  colnames(mtx.coverage) <- s.names
  colnames(mtx.m) <- s.names
  colnames(mtx.beta) <- s.names
  rownames(mtx.coverage) <- as.character(ids)
  rownames(mtx.m) <- as.character(ids)
  rownames(mtx.beta) <-as.character(ids)
  
  return (list("coverage"=mtx.coverage, "beta"=mtx.beta, "m"=mtx.m))
}

# get methylation data in a genomics range
getMethInGRange <- function(meth, gr) {
  meth.gr <- as(meth, "GRanges")
  hits <- !is.na(findOverlaps(meth.gr, gr, select = "first"))
  return (meth[hits, ])
}

##
plot_DMR_Gviz <- function(dmrs, dmrTSS, id, methDiff = myDiff.DSS) {
  # All DMR's associated with the gene of id
  dmr_id <- dmrTSS[dmrTSS$ID == id, ]
  # Get GRanges for the DMR's
  dmr_id.gr <- as(dmrs[dmr_id$target.row, ], "GRanges")
  strand(dmr_id.gr) <- "*"
  
  # Get exon structure of gene of id
  rna <- rownames(mapping[mapping$ID == id, ])
  exons <- gene.parts$exons[gene.parts$exons$name == rna, ]
  grlist <- GRangesList("exons" = exons[, 2], "dmr" = dmr_id.gr[, 1])
  gr <- unlist(grlist)
  start <- min(start(range(gr)))
  end <- max(end(range(gr)))
  chr <- as.character(seqnames(gr))[1]   # convert factor to string
  
  
  gtrack <- GenomeAxisTrack()
  #meth_id <- getMethInGRange(methDiff, range(gr))
  meth_id <- getMethInGRange(methDiff, range(dmr_id.gr))
  meth_id.gr <- as(meth_id, "GRanges")
  meth_id.dt <- DataTrack(range=meth_id.gr, data=meth_id.gr$meth.diff, legend = TRUE, start = start, end = end, type = "h", chromosome = chr, name = "Methyl Diff")
  exons.track <- AnnotationTrack(exons, name = "Exons", chromosome = chr, start = start, end = end, id = "exons")
  plotTracks(list(meth_id.dt, gtrack, exons.track), main = paste("Gene", id, chr), sizes=c(1, 0.5, 0.3), from = start, to = end)
}

##
# Plot DMR's associated with a gene of id
plotOneDMR <- function(dmrs, dmrTSS, id, methDiff = myDiff.DSS) {
  dmr_id <- dmrTSS[dmrTSS$ID == id, ]
  dmr_id.gr <- as(dmrs[dmr_id$target.row, ], "GRanges")
  gr <- range(dmr_id.gr)
  strand(gr) <- "*"
  methhits <- getMethInGRange(methDiff, gr)
  matplot(methhits$start, methhits$meth.diff, type="o", pch = 20,col = 'blue',
          xlab=paste("Genomic location", seqnames(gr), sep = " "), 
          ylab="Methylation difference" )
  abline(h=c(-10,0,10),lty=2)
}

plot_DMRs <- function(dmrs) {
  plot(dmrs$meth.diff, dmrs$logFC, type="p", cex = 1, pch = 20, col= "blue", xlab = "Methylation Difference", ylab = "Log Fold Change")
  text(dmrs$meth.diff, dmrs$logFC + 0.25, labels = dmrs$ID, cex = 0.5)
}

# Plot methylation data for a gene
plot_Methyl_Grange <- function(methData, grang) {
  # 
  methhits <- getMethInGRange(methData, grang)
  matplot(methhits$start, methhits$meth.diff, type="o", pch = 20,col = 'blue',
          xlab=paste("Genomic location", methhits$chr[1]), 
          ylab="Methylation difference" )
  abline(h=c(-10,0,10),lty=2)
}