#'@title GGI Matrix Plot
#'
#'@description \code{GGI.plot} is a graphic function that allows to visualize a
#'  GGI results matrix in a style similar to a heatmap. The plot possess an
#'  interactive functionality when p-values or genes names couldn't be drawn.
#'
#'@details This function draw the upper half of a GGI results matrix without its
#'  diagonal (as it holds no meaning). A gradient is created from 0 to 1 (by
#'  default from crimson to white) and the matrix cells are colored according to
#'  the corresponding p-value.
#'
#'  By  default, p-values and names are drawn to make matrix reading easier, but
#'  in case \code{GGI} matrix is large, p-values (and eventually genes' names as
#'  \code{GGI} grows bigger) are not drawn. In that case, the default behavior
#'  of the function is to start an interactive process where user can click on a
#'  cell of interest to open a tooltip displaying which genes are involved in
#'  selected interaction and the p-value of the interaction test. Tooltips can
#'  be closed if user clicks anywhere else than on a cell. This process stops
#'  when the user presses the escape button (or terminates the locator procedure
#'  in general) or when the user clicks on any place other than a cell when no
#'  tooltip window is open.
#'
#'  To improve plot clarity, user may set a threshold above which cells are
#'  colored with a distinct color. By default, threshold is set to 1 and no
#'  cell is colored differently (as values must be \emph{strictly} above the
#'  threshold).
#'
#'@param GGI Numeric symmetric matrix of size \eqn{G*G} where \eqn{G} is the
#'  number of genes involved in the GGI analysis. Matrix' dimensions may be
#'  named, if not generic names will be generated. \eqn{G} must be at least 3.
#'@param col String vector. Marker colors to be used for the gradient. The first
#'  element of the vector is the value for 0 and the last is for 1. If only one
#'  value is parsed, that color is used for 0 and white is automatically used
#'  for 1. Any value compatible with \code{\link{colorRampPalette}} function.
#'@param colbar.width A positive number describing the gradient bar width ratio.
#'  That number is used to keep the gradient bar's width steady as the size of
#'  \code{GGI} increases.
#'@param title A string used as the plot title. If left as NUUL, a generic name
#'  is generated.
#'@param hclust.order A boolean. Should a hierachical clustering procedure be
#'  performed on \code{GGI} to order the matrix ?
#'@param threshold A numeric between 0 and 1. All p-value strictly greater than
#'  that value are distinctly colored (See \code{NA.col}).
#'@param NA.col A string. The color to use when a p-value is strictly
#'  greater than \code{threshold}.
#'@param draw.pvals A boolean. Should p-values be plotted ? Disabled when
#'  \code{GGI}'s size exceeds \eqn{15*15}.
#'@param draw.names A boolean. Should genes' names be plotted on matrix margins
#'  ? Disabled when \code{GGI}'s size exceeds \eqn{25*25}.
#'@param interact A boolean. Should the plot be clickable ? (See Details for
#'  more information). Disabled when open R session is not interactive.
#'
#'@return Nothing is returned if the plot could be successfully plotted.
#'
#'  A warning is issued if \code{draw.names} and/or \code{draw.pvals} are
#'  manually set to TRUE when \code{GGI} is too big (respectively more than 25
#'  or 15 genes.)
#'
#'@seealso \code{\link{GGI}}
#'
#'@export
#'
#'@examples
#' ## Dataset is included in the package
#' ped <- system.file("extdata/example.ped", package="GGItest")
#' info <- system.file("extdata/example.info", package="GGItest")
#' posi <- system.file("extdata/example.txt", package="GGItest")
#'
#' dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\t")
#' dta <- imputeSnpMatrix(dta$snpX, genes.info = dta$genes.info)
#' resp <- system.file("extdata/response.txt", package="GGItest")
#'
#' # If a data frame is provided to GGI or one of the *.test function, only the
#' # column is checked and used.
#' Y  <- read.csv(resp, header=FALSE)
#'
#' ## By default the PCA-based method is used.
#' GGI <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info)
#'
#' GGI.plot(GGI, hclust.order=TRUE, threshold=0.05)
GGI.plot <- function(GGI, col=c("#D6604D", "#104E8B"), colbar.width=0.15,
                     title=NULL, hclust.order=FALSE, use.log=FALSE,
                     threshold=NULL, NA.col="#D3D3D3",
                     draw.pvals=(ncol(GGI) <= 15), draw.names=(ncol(GGI) <= 25),
                     interact=!(draw.pvals && draw.names)) {

  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.character(col)) {
    stop("col argument should be a character vector.")
  } else if (!is.numeric(colbar.width) || colbar.width < 0) {
    stop("colbar.width argument should be a positive numeric")
  } else if (!is.character(title) && !is.null(title)) {
    stop("title argument should be a string.")
  } else if (!is.logical(draw.pvals) | !is.logical(draw.names)) {
    stop("show.pvals & draw.names arguments should be logical.")
  } else if (!is.logical(hclust.order)) {
    stop("hclust.order argument should be logical.")
  } else if (!is.logical(use.log)) {
    stop("use.log argument should be logical")
  } else if (!is.logical(interact)) {
    stop("interact argument should be logical")
  } else if (!is.null(threshold) && is.numeric(threshold) && (threshold > 1 || threshold < 0)) {
    stop("threshold argument can not be a numeric greater than 1 or lesser than 0.")
  } else if (!is.null(threshold) && is.character(threshold) && threshold != "R") {
    stop("threshold argument can not be any other string than 'R'.")
  } else if (!is.character(NA.col)) {
    stop("NA.col argument should be a character.")
  }

  R.thresh <- c(0.001, 0.01, 0.05, 0.1)

  # If only one color is parsed, white is
  # used to complete the scale
  if (length(col) < 2) {
    col <- c(col, "#FFFFFF")
  }

  if (use.log){
    GGI <- -log10(GGI)
    diag(GGI) <- min(GGI[row(GGI) != col(GGI)])
    col <- rev(col)

    if (!is.null(threshold)) threshold <- -log10(threshold);

    R.thresh <- -log10(R.thresh)
  }

  col.FUN <- grDevices::colorRampPalette(col)

  # Names checking (generated if none)
  if (is.null(dimnames(GGI))){
    genes.names <- paste("Gene", 1:ncol(GGI), sep=".")
    dimnames(GGI) <- list(genes.names, genes.names)
  }

  # Clustering
  if (hclust.order) {
    GGI.clust <- hclust(as.dist(GGI))
    GGI <- GGI[GGI.clust$order, GGI.clust$order]
  }

  # Calculating plot size
  plot.setup(GGI, colbar.width, draw.names, threshold)

  # Draw color map
  rect.pos <- draw.matrix(GGI, col.FUN, threshold, NA.col, R.thresh, use.log)

  # Draw color legend bar
  leg <- draw.colbar(GGI, col.FUN, colbar.width, threshold, R.thresh, NA.col, use.log)

  if (!is.null(leg)) leg <- leg$rect

  # Draw genes names
  if (draw.names && ncol(GGI) <= 25) {
    draw.genes.names(dimnames(GGI), rect.pos)
  } else if (draw.names) {
    warning("GGI object is too big (26+ genes): genes names were not plotted.
            The use of the tooltip functionality is recommanded.")
  }

  # Draw p-values
  if (draw.pvals && ncol(GGI) <= 15) {
    draw.interp(GGI, rect.pos)
  } else if (draw.pvals) {
    warning("GGI object is too big (16+ genes): p-values were not plotted.
             The use of the tooltip functionality is recommanded.")
  }

  # Draw title
  if (is.null(title)) {
    title <- "Genes Interactions Matrix Plot"
  }
  title(main=title)

  # Activate tooltip functionality
  if (interact && interactive()) {
    writeLines("Click on a cell to get info on that cell.\nPress Esc. to leave.")

	  prime.plot <- recordPlot()
  	inter.tooltip  <- FALSE
  	keep.on    <- TRUE
  	while(keep.on) {
  	  coords <- locator(n=1)
  	  # If user forcefully stop the locator function then
  	  # break out of the loop.
  	  if (is.null(coords)) break
  	  # As the bottom left point is used to identify a square
  	  # on the plot, coordinates are floored.
  	  coords <- floor(as.numeric(coords))

  	  # Plot coordinates are converted back to matrix coordinates.
  	  coords <- c(row=nrow(GGI) - coords[2], col=coords[1])

  	  # Check if coordinates are conformant with GGI matrix
  	  coords.check <- try(GGI[coords['row'], coords['col']], silent=TRUE)
  	  if (class(coords.check) != 'try-error' && length(coords.check) == 1) {
  	    # Check if selected point is in upper triangle -diag excluded-
  	    # (onscreen part of the matrix).
  	    # It is the case when column index is strictly superior to
  	    # row index.
  	    if (coords[2] > coords[1]) {
  	      inter.tooltip <- TRUE
  	      clear.tooltip(inter.tooltip, prime.plot)
  	      draw.tooltip(coords, GGI, leg)
  	    } else {
  	      keep.on <- clear.tooltip(inter.tooltip, prime.plot)
  	      inter.tooltip <- !keep.on
  	    }
  	  } else {
  	    keep.on <- clear.tooltip(inter.tooltip, prime.plot)
  	    inter.tooltip <- !keep.on
  	  }
  	}
  }
}

# Function that computes the graphic window x and y extreme values.
# No real plotting happens in this function (empty window is opened).
plot.setup <- function(GGI, colbar.width, draw.names, threshold) {
  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.numeric(colbar.width)) {
    stop("colbar.width argument should be numeric.")
  } else if (!is.logical(draw.names)) {
    stop("draw.names argument should be TRUE or FALSE.")
  } else if (!is.null(threshold) && is.character(threshold) && threshold != "R") {
    stop("threshold argument can not be any other string than 'R'.")
  }

  # Widths and heights of elements are calculated
  if (is.null(threshold)) {
    colorLegend.height <- nrow(GGI)
    colorLegend.width  <- max(0.5, (ncol(GGI)/2)*colbar.width)
    matCol.padding <- colorLegend.width * 0.5
    colorLegend.space <- 1
  } else {
    colorLegend.height <- 0
    colorLegend.width <- 0.5
    matCol.padding <- colorLegend.width * 0.5
    colorLegend.space <- 0
  }

  plot.width  <- ncol(GGI) + colorLegend.space + colorLegend.width + matCol.padding
  plot.height <- nrow(GGI)

  plot(0, xlim=c(2, plot.width), ylim=c(1, plot.height), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n")

  # If names are to be plotted and exist, text padding is calculated
  if (draw.names & ncol(GGI) <= 25 & !is.null(colnames(GGI))){
    names.length <- strwidth(colnames(GGI))

    text.vpadding <- ceiling(max(sin(pi/4) * names.length[-ncol(GGI)])) + 0.25
    text.lpadding <- floor(min(seq(2, nrow(GGI)) - 0.25 - names.length[-1]))
    text.rpadding <- ceiling(max(cos(pi/4) * names.length[-ncol(GGI)]))
  } else {
    text.vpadding <- 2
    text.lpadding <- 0
    text.rpadding <- 0
  }

  if (is.null(threshold)) {
    colbar.text.padding <- ceiling(colorLegend.width*0.1 + strwidth("0.75"))
  } else {
    colbar.text.padding <- 0
  }

  xlim <- c(text.lpadding, plot.width + colbar.text.padding + text.rpadding)
  ylim <- c(1, plot.height + text.vpadding)
  plot(0, xlim=xlim, ylim=ylim, type="n",
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
}

# Function that draw the upper triangle of GGI matrix.
# Cells are colored according to corresponding p-values and
# the value of threshold.
# Invisibly return the coordinates of the bottom left point of each
# square. (to save some computing time later)
draw.matrix <- function(GGI, col.FUN, threshold, NA.col, R.thresh, use.log = FALSE) {
  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.function(col.FUN)) {
    stop("col.FUN argument should be a function resulting from colorRampPalette.")
  } else if (!is.logical(use.log)) {
    stop("use.log argument should be a logical.")
  } else if (!is.numeric(R.thresh) || length(R.thresh) != 4) {
    stop("R.thresh argument should be a numeric vector of length 4.")
  } else if (!is.character(NA.col)) {
    stop("NA.col argument should be a character.")
  }

  rect.data <- GGI[upper.tri(GGI)]

  # Assigning colors depending on display options
  if (is.null(threshold)){
    # If gradient is displayed then probs are turned into a percentage first
    if (use.log) {
      quantiles <- c(0, max(GGI, na.rm=TRUE))
    } else {
      quantiles <- c(0, 1)
    }

    rect.perc <- (rect.data - quantiles[1]) / (diff(quantiles))

  } else if (is.numeric(threshold)) {
    # If a threshold is used
    if (use.log) {
      rect.perc <- ifelse(rect.data >= threshold, 0, 1)
    } else {
      rect.perc <- ifelse(rect.data <= threshold, 0, 1)
    }
  } else if (is.character(threshold)) {
    if (use.log){
      rect.perc <- findInterval(rect.data, rev(R.thresh))

    } else {
      rect.perc <- findInterval(rect.data, R.thresh)
    }

      rect.perc <- rect.perc/4
  }

  rect.perc <- floor(rect.perc*200)
  rect.perc[rect.perc == 0] <- 1
  rect.col <- col.FUN(200)[rect.perc]

  # NA values are also disabled
  rect.col[which(is.na(rect.data))] <- NA.col

  rect.pos  <- which(upper.tri(GGI), arr.ind = TRUE)
  temp.X <- rect.pos[, 2]
  rect.pos[, 2] <- max(rect.pos[, 1]) - rect.pos[, 1] + 1
  rect.pos[, 1] <- temp.X

  rect(xleft = rect.pos[, 1],
       ybottom = rect.pos[, 2],
       xright  = rect.pos[, 1] + 1,
       ytop = rect.pos[, 2] + 1,
       col  = rect.col)

  invisible(rect.pos)
}

# Function that draws the gradient indicator.
# The gradient bar is sliced accross the height into a
# large number of smaller rectangles.
draw.colbar <- function(GGI, col.FUN, colbar.width, threshold, R.thresh, NA.col, use.log = FALSE) {
  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.function(col.FUN)) {
    stop("col.FUN argument should be a function resulting from colorRampPalette.")
  } else if (!is.logical(use.log)) {
    stop("use.log argument should be a logical.")
  } else if (!is.numeric(R.thresh) || length(R.thresh) != 4) {
    stop("R.thresh argument should be a numeric vector of length 4.")
  } else if (!is.character(NA.col)) {
    stop("NA.col argument should be a character.")
  }

  if (is.null(threshold)){
    colorLegend.height <- nrow(GGI)
    colorLegend.width  <- max(0.5, (ncol(GGI)/2)*colbar.width)
    matCol.padding <- colorLegend.width * 0.5

    NA.height <- 0.05 * colorLegend.height
    NA.padding <- 0.5 * matCol.padding
    colbar.start <- NA.height + NA.padding

    rect(xleft = rep(ncol(GGI) + 1 + matCol.padding, 200),
         ybottom = seq(1 + colbar.start, ncol(GGI), length=201)[-201],
         xright = rep(ncol(GGI) + 1 + matCol.padding + colorLegend.width, 200),
         ytop = seq(1 + colbar.start, ncol(GGI), length=201)[-1],
         col = col.FUN(200),
         border = NA)

    rect(xleft = ncol(GGI) + 1 + matCol.padding,
         ybottom = 1 + colbar.start,
         xright = ncol(GGI) + 1 + matCol.padding + colorLegend.width,
         ytop = ncol(GGI),
         col = NA,
         border = "black")

    if (use.log) {
      quantiles <- as.numeric(format(quantile(GGI, na.rm=TRUE, names=FALSE), digits=2))
      quantiles.pos <- seq(0, 1, 0.25)
    } else {
      quantiles <- c(0, 0.05, 0.25, 0.5, 0.75, 1)
      quantiles.pos <- quantiles
    }

    segments(x0 = rep(ncol(GGI) + 1 + matCol.padding + colorLegend.width, 6),
             y0 = (ncol(GGI) - 1 -colbar.start) * quantiles.pos + 1 + colbar.start,
             x1 = rep(ncol(GGI) + 1 + matCol.padding + colorLegend.width*1.1 , 6),
             y1 = (ncol(GGI) - 1 -colbar.start) * quantiles.pos + 1 + colbar.start)

    text(x = rep(ncol(GGI) + 1 + matCol.padding + colorLegend.width*1.1 , 6),
         y = (ncol(GGI) - 1 - colbar.start) * quantiles.pos + 1 + colbar.start,
         labels = quantiles,
         pos = 4)

    # NA legend
    rect(xleft = ncol(GGI) + 1 + matCol.padding,
         ybottom = 1,
         xright = ncol(GGI) + 1 + matCol.padding + colorLegend.width,
         ytop =  1 +NA.height,
         col = NA.col,
         border = "black")

    text(x = ncol(GGI) + 1 + matCol.padding + colorLegend.width*1.1,
         y = 1 + 0.5*NA.height,
         labels = "NA",
         pos = 4)

    return(NULL)

  } else if (is.numeric(threshold)) {
    sign <- c("<", ">")

    leg <- legend("bottomleft", c(paste(sign, round(threshold, 3)), 'NA'),
           fill = c(col.FUN(200)[c(1, 200)], NA.col)
           )
  } else if (is.character(threshold)) {
    if (!use.log) {
      legends <- c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "> 0.1")
    } else {
      legends <- c("< 1", "> 1", "> 1.3", "> 2", "> 3")
    }

    leg <- legend("bottomleft", c(legends, 'NA'),
           fill = c(col.FUN(200)[c(1, 50, 100, 150, 200)], NA.col)
    )
  }

  return(leg)
}

# Function that draws genes' names on the plot.
# As the diagonale of the matrix is not drawn, the first
# names is skipped vertically and the last horizontally.
draw.genes.names <- function(genes.names, rect.pos) {
  if(!is.list(genes.names) && length(genes.names) != 2) {
    stop("genes.names argument should be a list of length two.")
  } else if (!is.character(genes.names[[1]]) | !is.character(genes.names[[2]])) {
    stop("genes.names argument should be a list of character vectors.")
  } else if (length(genes.names[[1]]) != length(genes.names[[2]])) {
    stop("genes.names[[1]] & genes.names[[2]] should be of same length.")
  } else if (length(genes.names[[1]]) > 25) {
    stop("Can't handle more than 25 names.")
  } else if (!is.matrix(rect.pos) && !is.numeric(rect.pos[1, 1])) {
    stop("rect.pos argument should be a numeric matrix")
  }

  cex=sort((1 - 1/3:25), decreasing=TRUE)[length(genes.names[[1]]) - 2]

  # Horizontaly
  text(x = sort(unique(rect.pos[, 1])) - 0.25,
       y = sort(unique(rect.pos[, 2]), decreasing = TRUE) + 0.5,
       labels = genes.names[[1]][-length(genes.names[[2]])],
       pos = 2)

  # Verticaly
  text(x = sort(unique(rect.pos[, 1])) + 0.5*min(c(1, (1/length(genes.names[[1]]) ))),
       y = max(rect.pos[, 2]) + 1 + 0.25,
       labels = genes.names[[2]][-1],
       pos = 4,
       srt = 45)
}

# Function that draws the interaction p-values on the matrix.
# Multiple cex are tested so that it is ensured that p-values
# fit inside the squares and are still big enough.
draw.interp <- function(GGI, rect.pos){
  if(!is.matrix(GGI) && !is.numeric(GGI[1, 1])) {
    stop("GGI argument should be a numeric matrix.")
  } else if (ncol(GGI) != nrow(GGI)) {
    stop("GGI argument should a symmetric matrix.")
  } else if (ncol(GGI) < 3) {
    stop("At least 3 genes must be provided.")
  } else if (!is.matrix(rect.pos) && !is.numeric(rect.pos[1, 1])) {
    stop("rect.pos argument should be a numeric matrix")
  }

  rect.data <- GGI[upper.tri(GGI)]
  rect.data <- format(rect.data, digits=2, scientific=TRUE)

  for (i in seq(1, 0, length=30)) {
    cex <- i
    pval.width <- strwidth(rect.data, cex=cex)
    if (max(pval.width) < 0.9) { break }
  }

  x <- rect.pos[, 1] + 0.5
  y <- rect.pos[, 2] + 0.5
  text(x, y, labels=rect.data, cex=cex)
}

# Function that draws the tooltip windows
# A white black-bordered box is first created
# and text is plotted on top of it.
draw.tooltip <- function(coords, GGI, legend.box) {
  if (is.null(legend.box)) {
    bottomleft <- par('usr')[c(1, 3)]
  } else  {
    bottomleft <- c(legend.box$left + legend.box$w, legend.box$top - legend.box$h)
  }

  tooltip.str <- paste0('Interaction: ',
                       rownames(GGI)[coords[1]], ':',
                       colnames(GGI)[coords[2]],
                       '\np-val: ',
                       format(GGI[coords[1], coords[2]],
                              digits=4))

  rect(xleft = bottomleft[1],
       ybottom = bottomleft[2],
       xright = bottomleft[1] + strwidth(tooltip.str)*1.1,
       ytop = bottomleft[2] + strheight(tooltip.str)*1.5,
       col = "white")
  text(x = bottomleft[1] + strwidth(tooltip.str)*0.05,
       y = mean(c(bottomleft[2], bottomleft[2] + strheight(tooltip.str)*1.5)),
       labels = tooltip.str,
       pos = 4, offset=0)

}

# Function that handles tooltip clearing and tooltip
# procedure ending.
clear.tooltip <- function(inter.tip, prime.plot) {
  # If not in upper triangle then tooltip if cleared
  if (inter.tip) {
    replayPlot(prime.plot)
    return(TRUE)
  } else {
    # If tooltip is already cleared then interaction with
    # user is ceased.
    return(FALSE)
  }
}
