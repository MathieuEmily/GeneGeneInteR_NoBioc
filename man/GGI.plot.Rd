\name{GGI.plot}
\alias{GGI.plot}
\title{Plot a Gene-Gene Interaction Study with an heatmap-like visualization}
\usage{
GGI.plot(GGI, col = c("#D6604D", "#104E8B"), colbar.width = 0.15,
  title = NULL, hclust.order = FALSE, use.log = FALSE, threshold = NULL,
  NA.col = "#D3D3D3", draw.pvals = (ncol(GGI) <= 15),
  draw.names = (ncol(GGI) <= 25), interact = !(draw.pvals && draw.names))
}
\arguments{
\item{GGI}{A \code{\link{matrix}} or a \code{\link{data.frame}} of pValues containing the results of a set of genes interaction analysis. \code{GGI} must be a squared matrix of size \eqn{G*G} where \eqn{G} is the number of genes considered. \code{GGI} can be produced as the output of the function \code{\link{GGI}}. Matrix' dimensions may be named, if not generic names will be generated. \eqn{G} must be at least 3.}

\item{col}{String vector. Marker colors to be used for the gradient. The first
element of the vector is the value for 0 and the last is for 1. If only one
value is parsed, that color is used for 0 and white is automatically used
for 1. Any value compatible with \code{\link{colorRampPalette}} function.}

\item{colbar.width}{A positive number describing the gradient bar width ratio.
That number is used to keep the gradient bar's width steady as the size of
\code{GGI} increases.}

\item{title}{A string used as the plot title. If left as NULL, a generic name
is generated.}

\item{hclust.order}{A boolean. Should a hierachical clustering procedure be
performed on \code{GGI} to order the matrix ?}

\item{threshold}{A numeric between 0 and 1. All p-value strictly greater than
that \code{threshold} are distinctly colored (See \code{NA.col}).}

\item{NA.col}{A string. The color used when a p-value is strictly
greater than \code{threshold}.}

\item{draw.pvals}{A boolean. Should p-values be plotted ? Disabled when
the size of \code{GGI} exceeds \eqn{15*15}.}

\item{draw.names}{A boolean. Should genes' names be plotted on matrix margins
? Disabled when the size of \code{GGI} exceeds \eqn{25*25}.}

\item{interact}{A boolean. Should the plot be clickable ? (See Details for
more information). Disabled when open R session is not interactive.}
}
\description{
\code{GGI.plot} is a graphical function that allow the heatmap-like visualization of a Gene-Gene Interaction study based on a set of genes.
 The plot possess an interactive functionality when p-values or genes names couldn't be drawn.
}
\details{
This function draw the upper half of a Gene-Gene Interaction results matrix without its
 diagonal (as it holds no meaning). A gradient is created from 0 to 1 (by
 default from crimson to white) and the matrix cells are colored according to
 the corresponding p-value.

 By  default, p-values and names are drawn to make matrix reading easier, but
 in case parameter \code{GGI} is large, p-values (and eventually gene names as
 \code{GGI} grows bigger) are not drawn. In that case, the default behavior
 of the function is to start an interactive process where user can click on a
 cell of interest to open a tooltip displaying which genes are involved in
 selected interaction and the p-value of the interaction test. Tooltips can
 be closed if user clicks anywhere else than on a cell. This process stops
 when the user presses the escape button (or terminates the locator procedure
 in general) or when the user clicks on any place other than a cell when no
 tooltip window is open.

 To improve plot clarity, user may set a threshold above which cells are
 colored with a distinct color. By default, threshold is set to 1 and no
 cell is colored differently (as values must be \emph{strictly} above the
 threshold).
}
\examples{
## Dataset is included in the package
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")

dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\\t")
dta <- imputeSnpMatrix(dta$snpX, genes.info = dta$genes.info)
resp <- system.file("extdata/response.txt", package="GeneGeneInteR")

# If a data frame is provided to GGI or one of the *.test function, only the
# column is checked and used.
Y  <- read.csv(resp, header=FALSE)

## By default the PCA-based method is used.
GGI.res <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info)

GGI.plot(GGI.res, hclust.order=TRUE, threshold=0.05)
}
\seealso{
\code{\link{GGI}}
}

