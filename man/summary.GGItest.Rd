\name{summary.GGItest}
\alias{summary.GGItest}
\title{Printing summaries class \code{"GGItest"} objects}
\usage{
	## S3 method for class 'GGItest'
summary.GGItest(x)
}
\arguments{
\item{x}{\code{GGItest}. Given as output of \code{\link{CCA.test}}.}

}
\description{
	Printing objects of class \code{"GGItest"}
}
\details{
\code{print.GGItest} provides a customized summary of a \code{GGItest} object.
}
\examples{

## Not run:
data(gene.pair)
summary(PCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,threshold=0.7,method="Std"))

## End (Not run)
}

