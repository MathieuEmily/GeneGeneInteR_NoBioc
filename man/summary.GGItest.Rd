\name{summary.GGItest}
\alias{summary.GGItest}
\title{Printing summaries class \code{"GGItest"} objects}
\usage{
	## S3 method for class 'GGItest'
\method{summary}{GGItest}(object, \dots)
}
\arguments{
\item{object}{\code{GGItest}. Given as output of \code{\link{CCA.test}}.}
\item{\dots}{further arguments passed to or from other methods}
}
\description{
	Printing summaries for objects of class \code{"GGItest"}
}
\details{
\code{summary.GGItest} provides a customized summary of a \code{GGItest} object.
}
\examples{

## Not run:
## Dataset is included in the package
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")

## Importation of the genotypes
dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\\t")
## Filtering of the data: SNPs with MAF < 0.05 or p.value for HWE < 1e-3 are removed. 
# No filtering is applied regarding missing data (call.rate=0).
dta <- snpMatrixScour(snpX=dta$snpX,genes.info=dta$genes.info,min.maf=0.05,min.eq=1e-3,call.rate=0)
## Imputation of the missing genotypes
dta <- imputeSnpMatrix(dta$snpX, genes.info = dta$genes.info)

## Importation of the phenotype
resp <- system.file("extdata/response.txt", package="GeneGeneInteR")
Y  <- read.csv(resp, header=FALSE)

## plot of the interaction between the 17 genes with the CLD method -- can take a few minutes
\dontrun{
GGI.res <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info,method="CLD")
GGI.plot(GGI.res,threshold=0.05)
}

## plot of the interaction between 12 among the 17 genes with the default PCA method 
## Selection of 12 genes among 17
dta <- select.snps(dta$snpX, dta$genes.info, c("bub3","CDSN","Gc","GLRX","PADI1","PADI2","PADI4",
                   "PADI6","PRKD3","PSORS1C1","SERPINA1","SORBS1"))
GGI.res <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info)

## Summary of the results with default values
summary(GGI.res)

## End (Not run)
}

