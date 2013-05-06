#' Functions for analyzing and visualizing viral integration site data.
#'
#' hiAnalyzer contains set of functions which allow users to analyze and visualize viral integration site data. Most of the data used by this package is an output of either hiAnnotator or hiReadsProcessor package.
#'
#' @import ggbio scales gridExtra grid seqLogo RColorBrewer reshape dataframe
#' @docType package
#' @name hiAnalyzer
NULL

#### genome plot ####
#' @export
plot.CircleIdeogram <- function(sites, separateByCol = "Alias", freeze = "hg18", pointSize = 0.5, 
                                yAxisScatter = NULL, title = "", clusterSitesBin = NULL, adjustTracks = TRUE, dbcon = NULL) {
  stopifnot(nrow(sites) > 0)
  require(hiAnnotator)
  
  dbconDisconnect <- FALSE
  if (is.null(dbcon)) {
    dbconDisconnect <- TRUE
    dbcon <- connectToIntSitesDB()
  }
  
  # get tracks #
  sites.gr <- as(makeRangedData(subset(sites, type == "insertion"), soloStart = TRUE, 
                                chromCol = "Chr"), "GRanges")
  genes <- as(makeRangedData(unique(dbGetQuery(dbcon, paste("select * from ", freeze, 
                                                            ".refflat", sep = ""))[, 3:6])), "GRanges")
  freeze.gr <- as(makeRangedData(unique(dbGetQuery(dbcon, paste("select chrom, 1 as start, size as end ,'*' as strand from ", 
                                                                freeze, ".chrominfo", sep = "")))), "GRanges")
  seqlengths(freeze.gr) <- end(freeze.gr)
  seqorder <- sapply(strsplit(gsub("chr", "", seqlevels(freeze.gr)), "\\_"), "[[", 
                     1)
  freeze.gr <- freeze.gr[order(as.numeric(seqorder))]
  seqlevels(freeze.gr) <- seqlevels(freeze.gr)[order(as.numeric(seqorder))]
  
  if (dbconDisconnect) {
    dbDisconnect(dbcon)
  }
  
  # reorder chromosome names, add sizes, and do basic formating if any #
  seqlevels(sites.gr) <- seqlevels(freeze.gr)[seqlevels(freeze.gr) %in% seqlevels(sites.gr)]
  seqlengths(sites.gr) <- seqlengths(freeze.gr)[seqlevels(freeze.gr) %in% seqlevels(sites.gr)]
  seqlevels(genes) <- seqlevels(freeze.gr)[seqlevels(freeze.gr) %in% seqlevels(genes)]
  seqlengths(genes) <- seqlengths(freeze.gr)[seqlevels(freeze.gr) %in% seqlevels(genes)]
  freeze.gr <- keepSeqlevels(freeze.gr, seqlevels(sites.gr))
  genes <- keepSeqlevels(genes, seqlevels(sites.gr))
  
  newnames <- structure(gsub("chr", "", seqlevels(freeze.gr)), names = seqlevels(freeze.gr))
  freeze.gr <- renameSeqlevels(freeze.gr, newnames)
  sites.gr <- renameSeqlevels(sites.gr, newnames)
  genes <- renameSeqlevels(genes, newnames)
  
  values(freeze.gr)$Chr <- as.character(seqnames(freeze.gr))
  
  # cluster genes to reduce amount of data getting plotted #
  genes <- genes[order(genes)]
  strand(genes) <- "*"
  genes.reduced <- reduce(genes, min.gapwidth = 1e+05)
  values(genes.reduced)$counts <- countOverlaps(genes.reduced, genes, ignore.strand = T)
  values(genes.reduced)$logCount <- log(values(genes.reduced)$counts)
  rm(genes)
  
  # cluster sites if defined #
  asBars <- FALSE
  if (!is.null(clusterSitesBin)) {
    asBars <- TRUE
    sites.list <- sites.gr
    strand(sites.list) <- "*"
    values(sites.list) <- NULL
    sites.list <- split(sites.list, values(sites.gr)[, separateByCol])
    sites.reduced <- GRanges()
    for (x in names(sites.list)) {
      reduced <- reduce(sites.list[[x]], min.gapwidth = clusterSitesBin)
      values(reduced)[, separateByCol] <- x
      values(reduced)$counts <- countOverlaps(reduced, sites.list[[x]], ignore.strand = T)
      sites.reduced <- c(sites.reduced, reduced)
    }
    sites.gr <- sites.reduced
    values(sites.gr)$logCount <- log(values(sites.gr)$counts)
    rm("sites.list", "sites.reduced", "reduced")
  }
  
  ## make the plot ##
  tracks <- sort(unique(values(sites.gr)[, separateByCol]))
  dataset <- structure(c("#FF0000", getColors(length(tracks))), names = c("Genes", tracks))
  
  if (adjustTracks) {
    track.radius <- structure(seq(23, pmax(23/length(tracks), 2), length = length(tracks)), 
                              names = tracks)
  } else {
    rads <- seq(23, 2, by = -4.5)
    if (length(rads) < length(tracks)) {
      stop("No enough space available to plot all the tracks. Try using adjustTracks=TRUE")
    }
    track.radius <- structure(rads[1:length(tracks)], names = tracks)
  }
  
  p <- ggplot() + layout_circle(freeze.gr, geom = "ideo", fill = "gray70", radius = 34, 
                                trackWidth = 2)
  p <- p + layout_circle(freeze.gr, geom = "scale", size = 2, radius = 36, trackWidth = 2)
  p <- p + layout_circle(freeze.gr, geom = "text", aes(label = Chr), vjust = 0, radius = 38, 
                         trackWidth = 7)
  p <- p + layout_circle(genes.reduced, aes(y = logCount), geom = "bar", fill = "#FF0000", 
                         color = "red", radius = 28, trackWidth = 6, size = 0.2)
  
  # make a track per separateByCol from sites.gr #
  for (f in tracks) {
    dat <- subset(sites.gr, values(sites.gr)[, separateByCol] == f)
    if (asBars) {
      p <- p + layout_circle(dat, aes(y = logCount), geom = "bar", fill = dataset[[f]], 
                             color = dataset[[f]], radius = track.radius[[f]], trackWidth = 6, size = 0.2)
    } else {
      if (is.null(yAxisScatter)) {
        values(dat)$YaxisScatter <- rnorm(length(dat), 5)
      } else {
        values(dat)$YaxisScatter <- values(dat)[, yAxisScatter]
      }
      p <- p + layout_circle(dat, aes(y = YaxisScatter), size = pointSize, geom = "point", 
                             color = dataset[[f]], radius = track.radius[[f]])
    }
  }
  p <- p + theme(title = title)
  
  test <- qplot(data = data.frame(Track = factor(names(dataset), levels = names(dataset))), 
                1, Track, fill = Track, geom = "tile") + scale_fill_manual(values = dataset)
  tmp <- ggplot_gtable(ggplot_build(test))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  
  grid.arrange(p, legend, widths = unit.c(unit(1, "npc") - legend$width, legend$width), 
               nrow = 1)
}

#### weblogo plot ####
#' @export
plot.siteLogo <- function(seqs, samples, traditional = FALSE, ZeroCentered = TRUE, themes = NULL) {
  seqs <- split(as.character(seqs), as.character(samples))
  
  dat <- sapply(seqs, function(x) consensusMatrix(DNAStringSet(x), as.prob = T)[1:4, 
                                                                                ], simplify = F)
  if (ZeroCentered) {
    dat <- lapply(dat, function(x) {
      newones <- seq(-(ncol(x)/2), (ncol(x)/2), by = 1)
      colnames(x) <- newones[newones != 0]
      x
    })
  }
  
  if (traditional) {
    sapply(names(dat), function(x) {
      quartz(x)
      seqLogo(makePWM(dat[[x]]))
    })
  } else {
    melted <- do.call(rbind, lapply(names(dat), function(x) {
      cbind(samples = x, melt(cbind(as.data.frame(dat[[x]]), base = rownames(dat[[x]])), 
                              id = "base"))
    }))
    melted$variable <- gsub("V", "", melted$variable)
    melted$variable <- factor(melted$variable, levels = as.character(sort(unique(as.numeric(melted$variable)))))
    p <- qplot(data = melted, variable, value, geom = "line", colour = base, group = base, 
               xlab = "Position") + scale_y_continuous("Percent of Sequences", label = percent) + 
      facet_wrap(~samples, ncol = 1)
    
    if (!is.null(themes)) {
      p + get(themes)()
    } else {
      p
    }
  }
}

#### Distance to TSS plot ####
#' @export
plot.distToFeature <- function(dat, brks, byCol, annotCol, ratioMRC = FALSE, discreteBins = TRUE, 
                               themes = NULL, printPlotData = FALSE, bw = FALSE) {
  counts.df <- as.data.frame(xtabs(~type + get(byCol), dat, drop = T))
  names(counts.df)[2:3] <- c(byCol, "Total")
  plot.frame <- as.data.frame(xtabs(~cut(dat[, annotCol], brks, include.lowest = T, 
                                         dig.lab = 5) + get(byCol) + type, dat, drop = T))
  names(plot.frame)[1:2] <- c("Distance", byCol)
  plot.frame <- merge(plot.frame, counts.df)
  plot.frame$Percent <- with(plot.frame, Freq/Total)
  plot.frame$Dist.To.TSS <- as.numeric(sub(".+,(.+)]", "\\1", plot.frame$Distance))
  if (ratioMRC) {
    plot.frame <- arrange(plot.frame, get(byCol), type, Dist.To.TSS)
    ratiosCalc <- ddply(plot.frame, ~get(byCol) + annot + Dist.To.TSS, summarise, 
                        RatioMRC = Percent[tolower(type) == "insertion"]/Percent[tolower(type) %in% 
                                                                                   c("match", "mrcs")])
    res <- droplevels(subset(merge(plot.frame, ratiosCalc, all.x = T), type == "insertion"))
    if (!discreteBins) {
      p <- qplot(data = res, Dist.To.TSS, RatioMRC, geom = "bar", stat = "identity", 
                 position = "dodge", fill = get(byCol)) + scale_x_continuous("Distance", 
                                                                             expand = c(0, 0))
    } else {
      p <- qplot(data = res, Distance, RatioMRC, geom = "bar", stat = "identity", 
                 position = "dodge", fill = get(byCol)) + scale_x_discrete("Distance", 
                                                                           expand = c(0, 0))
    }
    p <- p + scale_y_continuous("Ratio Insertions/Matches", expand = c(0, 0))
  } else {
    if (!discreteBins) {
      p <- qplot(data = plot.frame, Dist.To.TSS, Percent, geom = "bar", stat = "identity", 
                 position = "dodge", fill = get(byCol)) + scale_x_continuous("Distance", 
                                                                             expand = c(0, 0))
    } else {
      p <- qplot(data = plot.frame, Distance, Percent, geom = "bar", stat = "identity", 
                 position = "dodge", fill = get(byCol)) + scale_x_discrete("Distance", 
                                                                           expand = c(0, 0))
    }
    p <- p + scale_y_continuous("Percent of Sites", label = percent, expand = c(0, 
                                                                                0)) + facet_wrap(~type, ncol = 1)
  }
  
  if (!is.null(themes)) {
    p <- p + get(themes)() + theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                              vjust = 0.5))
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  if (length(as.character(unique(dat[, byCol]))) > 9) {
    dat.colors <- getColors(length(as.character(unique(dat[, byCol]))))
    p <- p + scale_fill_manual(name = byCol, values = dat.colors)
  } else {
    p <- p + scale_fill_discrete(name = byCol)
  }
  
  if (bw) {
    p <- p + scale_fill_grey(name = byCol, start = 0.3, end = 0.7)
  }
  
  if (printPlotData) {
    plot.frame$byCol <- plot.frame[, byCol]
    print(cast(data = plot.frame, Distance ~ byCol, value = "Percent", margins = TRUE, 
               fun = sum))
  }
  
  p
}

#### Number/Percent of sites in feature ####
#' @export
plot.inFeature <- function(dat, byCol, annotCol, ratioMRC = FALSE, themes = NULL, 
                           printPlotData = FALSE, bw = FALSE) {
  dat[,annotCol] <- as.logical(dat[,annotCol])
  if(any(is.na(dat[,annotCol]))) {
    stop("Values in '", annotCol, "'' column are not boolean or only TRUE/FALSE")
  }
  
  counts.df <- count(dat, c("type", byCol))
  names(counts.df)[3] <- "Total"
  plot.frame <- count(dat, c("type", byCol, annotCol))
  plot.frame <- merge(plot.frame, counts.df)
  plot.frame$Percent <- with(plot.frame, freq/Total)
  
  if (ratioMRC) {
    plot.frame <- arrange(plot.frame, get(byCol), type, get(annotCol))
    ratiosCalc <- ddply(plot.frame, .(get(byCol),get(annotCol)), summarise, 
                        RatioMRC = Percent[tolower(type) == "insertion"]/
                          Percent[tolower(type) %in% c("match", "mrcs")])
    names(ratiosCalc)[1:2] <- c(byCol, annotCol)
    res <- droplevels(subset(merge(plot.frame, ratiosCalc, all.x = T), 
                             type == "insertion"))
    p <- qplot(data = res, get(annotCol), RatioMRC, geom = "bar", 
               stat = "identity", position = "dodge", fill = get(byCol)) +
      scale_y_continuous("Ratio Insertions/Matches", expand = c(0, 0)) +
      geom_hline(y=1)
      
  } else {
    p <- qplot(data = plot.frame, get(annotCol), Percent, geom = "bar", 
               stat = "identity", position = "dodge", fill = get(byCol)) + 
      scale_y_continuous("Percent of Sites", label = percent, 
                         expand = c(0, 0)) + facet_wrap(~type, ncol = 1)
  }
  
  p <- p + scale_x_discrete("In Gene", expand = c(0, 0))
  
  if (!is.null(themes)) {
    p <- p + get(themes)()      
  }
  
  if (length(as.character(unique(dat[, byCol]))) > 9) {
    dat.colors <- getColors(length(as.character(unique(dat[, byCol]))))
    p <- p + scale_fill_manual(name = byCol, values = dat.colors)
  } else {
    p <- p + scale_fill_discrete(name = byCol)
  }
  
  if (bw) {
    p <- p + scale_fill_grey(name = byCol, start = 0.3, end = 0.7)
  }
  
  if (printPlotData) {
    plot.frame$byCol <- plot.frame[, byCol]
    plot.frame$annotCol <- plot.frame[, annotCol]
    print(cast(data = plot.frame, annotCol ~ byCol, value = "Percent",
               margins = TRUE, 
               fun = sum))
  }
  
  p
}

### collection of brewer colors ###
#' @export
getColors <- function(n) {
  allCols <- unique(c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
                      brewer.pal(8, "Set2"), brewer.pal(12, "Paired"), brewer.pal(10, "Spectral")))
  if (n > length(allCols)) {
    stop("Not enough discrete colours to choose from! Max limit is ", length(allCols))
  }
  return(as.character(allCols[1:n]))
} 
