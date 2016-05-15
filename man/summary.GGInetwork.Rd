\name{summary.GGInetwork}
\alias{summary.GGInetwork}
\title{Printing summaries class \code{"GGInetwork"} objects}
\usage{
	## S3 method for class 'GGInetwork'
\method{summary}{GGInetwork}(object, \dots)
}
\arguments{
\item{object}{\code{GGInetwork}. Given as output of \code{\link{GGI}}.}
\item{\dots}{further arguments passed to or from other methods}
}
\description{
	Printing summaries for objects of class \code{"GGInetwork"}
}
\details{
\code{summary.GGItest} provides a customized summary of a \code{GGInetwork} object.
}
\examples{

## Not run:
data(gene.pair)
summary(PCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,threshold=0.7,method="Std"))

## End (Not run)
}

