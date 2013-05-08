#' Functions for analyzing and visualizing viral integration site or genomic data.
#'
#' hiAnalyzer contains set of functions which allow users to analyze and visualize viral integration site or genomic data. Most of the data used by this package is an output of either hiAnnotator or hiReadsProcessor package.
#'
#' @import ggbio scales gridExtra grid seqLogo RColorBrewer reshape hiAnnotator
#' @docType package
#' @name hiAnalyzer
NULL

#' Create a Circos styled genomic plot/ideogram.
#'
#' @param sites a GRanges object containing point genomic data with populated seqinfo slot. See \code{\link{makeGRanges}} for an easy way to do this!
#' @param separateByCol a column name in sites by which additional tracks are generated. Default is "group".
#' @param tracks a named list of GRanges object which serves as the target/reference/scaffold for data in sites object. Default is NULL. Example: list("genes"=genes.rd, "CpG"=cpgs.rd)
#' @param clusterSitesBin numeric bin size used to cluster data in sites. Default is NULL
#' @param clusterTracksBin numeric bin size used to cluster data in each object of tracks. Default is NULL
#' @param pointSize numeric value for the size of a point to denote a row in sites object. Default is 0.5. Not used if clusterSitesBin is defined, which in turn makes a bar plot.
#' @param title Title of the plot. Default is blank.
#'
#' @seealso \code{\link{plot.siteLogo}}, \code{\link{plot.distToFeature}}, 
#' \code{\link{plot.inFeature}}, \code{\link{getColors}}, \code{\link{makeGRanges}}
#'
#' @export
#'
#' @examples
plot.CircleIdeogram <- function(sites, separateByCol = "group", tracks = NULL, 
                                clusterSitesBin = NULL, clusterTracksBin = NULL,
                                pointSize = 0.5, title = "") {
  stopifnot(length(sites) > 0 | is(sites,"GRanges"))
  stopifnot(length(tracks) > 0 | all(sapply(tracks, is, class2="GRanges")))
  
  # cluster sites & tracks to reduce amount of data getting plotted #
  if(!is.null(clusterTracksBin)) {
    tracks <- sapply(tracks, 
                     function(x) {
                       x <- sort(x)
                       strand(x) <- "*"
                       x.reduced <- reduce(x, min.gapwidth = clusterTracksBin)
                       mcols(x.reduced)$counts <- countOverlaps(x.reduced, x, 
                                                                 ignore.strand = T)
                       mcols(x.reduced)$logCount <- log(mcols(x.reduced)$counts)
                       x.reduced
                     }, simplify=F)    
  }
  
  asBars <- FALSE
  if (!is.null(clusterSitesBin)) {
    asBars <- TRUE
    sites <- sapply(split(sites, mcols(sites)[[separateByCol]]), 
                    function(x) {
                      x <- sort(x)
                      strand(x) <- "*"
                      mcols(x) <- NULL
                      x.reduced <- reduce(x, min.gapwidth = clusterSitesBin)
                      mcols(x.reduced)$counts <- countOverlaps(x.reduced, x, 
                                                                ignore.strand = T)
                      mcols(x.reduced)$logCount <- log(mcols(x.reduced)$counts)
                      x.reduced
                    }, simplify=F) 
  }
  
  tracks <- c(tracks,sites)
  
  ## pick the object with the most seqnames to serve as the scaffold ##
  freeze.gr <- 
  
  ## make the plot scaffold ##
  dataset <- structure(getColors(length(tracks)), names = names(tracks))
  track.radius <- structure(seq(28, pmax(28/length(tracks), 2), length = length(tracks)), 
                            names = names(tracks))
  
  p <- ggplot()
  p <- p + layout_circle(freeze.gr, geom = "ideo", fill = "gray70", 
                         radius = 34, trackWidth = 2)
  p <- p + layout_circle(freeze.gr, geom = "scale", size = 2, 
                         radius = 36, trackWidth = 2)
  p <- p + layout_circle(freeze.gr, geom = "text", aes(label = Chr), vjust = 0, 
                         radius = 38, trackWidth = 7)
  
  # add tracks #
  for (f in tracks) {
    dat <- subset(sites.gr, mcols(sites.gr)[[separateByCol]] == f)
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
  
  test <- qplot(data = data.frame(Track = factor(names(dataset), 
                                                 levels = names(dataset))), 
                1, Track, fill = Track, geom = "tile") + 
    scale_fill_manual(values = dataset)
  tmp <- ggplot_gtable(ggplot_build(test))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  
  grid.arrange(p, legend, widths = unit.c(unit(1, "npc") - legend$width, legend$width), 
               nrow = 1)
}

#### weblogo plot ####
#' @export
plot.siteLogo <- function(seqs, samples, traditional = FALSE, ZeroCentered = TRUE, 
                          themes = NULL) {
  seqs <- split(as.character(seqs), as.character(samples))
  
  dat <- sapply(seqs, 
                function(x) consensusMatrix(DNAStringSet(x), as.prob = T)[1:4,], 
                simplify = F)
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
      cbind(samples = x, melt(cbind(as.data.frame(dat[[x]]), 
                                    base = rownames(dat[[x]])), 
                              id = "base"))
    }))
    melted$variable <- gsub("V", "", melted$variable)
    melted$variable <- factor(melted$variable, 
                              levels = as.character(sort(unique(as.numeric(melted$variable)))))
    p <- qplot(data = melted, variable, value, geom = "line", colour = base, 
               group = base, xlab = "Position") + 
      scale_y_continuous("Percent of Sequences", label = percent) + 
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
plot.distToFeature <- function(dat, brks, byCol, annotCol, ratioMRC = FALSE, 
                               discreteBins = TRUE, themes = NULL, 
                               printPlotData = FALSE, bw = FALSE) {
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
