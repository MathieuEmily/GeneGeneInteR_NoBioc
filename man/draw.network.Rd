\name{draw.network}
\alias{draw.network}
\title{Draw a Gene-Gene Interaction Network}
\usage{
draw.network(GGI, genes = 1:ncol(GGI),
  threshold = 0.05, plot.nointer = TRUE)
}
\arguments{
\item{GGI}{A \code{\link{matrix}} or a \code{\link{data.frame}} of pValues containing the results of a set of genes interaction analysis. \code{GGI} must be a squared matrix of size \eqn{G*G} where \eqn{G} is the number of genes considered. \code{GGI} can be produced as the output of the function \code{\link{GGI}}.}

\item{genes}{Numeric vector allowing a selection of the genes that will be included in the relations. Default is set to all genes.}

\item{threshold}{Numeric comprised in [0,1], defining the pValues that will be
considered as respresentative of an interaction. Default is set to 0.05.}

\item{plot.nointer}{A boolean. Set TRUE if the genes with no interaction
should be plotted.}
}
\description{
\code{\link{draw.network}} plot a graph representing the significative
interactions between genes of a Gene-Gene Interaction study.
}
\details{
The function takes as input the output of \code{\link{GGI}}, a matrix of pvalues
representing the level of interaction between genes. \code{genes} can be used
to select the genes you want to plot on the interaction network. It is either
a character vector giving genes' name or a numeric vector giving their column
number.
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
pVal <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info)

## By default, all genes are selected, threshold is 0.05.
draw.network(pVal)

## By default, genes with no interaction are plotted, use plot.nointer to change this option.
draw.network(pVal, threshold=0.005)
draw.network(pVal, threshold=0.005, plot.nointer=F)
}
\seealso{
\code{\link{draw.network}}, \code{\link{GGI}}
}

