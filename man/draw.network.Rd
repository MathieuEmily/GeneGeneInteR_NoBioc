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

## Importation of the genotypes
dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\\t")
## Filtering of the data: SNPs with MAF < 0.05 or p.value for HWE < 1e-3 are removed. No filtering is applied regarding missing data (NA.rate=1).
dta <- snpMatrixScour(snpX=dta$snpX,genes.info=dta$genes.info,min.maf=0.05,min.eq=1e-3,NA.rate=1)
## Imputation of the missing genotypes
dta <- imputeSnpMatrix(dta$snpX, genes.info = dta$genes.info)

## Importation of the phenotype
resp <- system.file("extdata/response.txt", package="GeneGeneInteR")
Y  <- read.csv(resp, header=FALSE)

## Network estimation of 17 genes with the CLD method -- can take a few minutes
\dontrun{
GGI.test <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info,method="CLD")
draw.network(pVal,threshold=0.05)
}

## Network estimation of 12 among the 17 genes with the default PCA method 
## Selection of 12 genes among 17
dta <- select.snps(dta$snpX, dta$genes.info, c("bub3","CDSN","Gc","GLRX","PADI1","PADI2","PADI4","PADI6","PRKD3","PSORS1C1","SERPINA1","SORBS1"))
GGI.test <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info)

## Example of network with default threshold 0.05
draw.network(GGI.test)

## Example of network with threshold 0.01 where genes with no interaction are not plotted (plot.nointer=FALSE)
draw.network(GGI.test,threshold=0.1,plot.nointer=FALSE)

}
\seealso{
\code{\link{draw.network}}, \code{\link{GGI}}
}

