\name{GBIGM.test}
\alias{GBIGM.test}
\title{GBIGM (Gene-based information gain method) for GGI analysis.}
\usage{
GBIGM.test(Y, G1, G2, n.perm = 1000)
}
\arguments{
\item{Y}{numeric or factor vector with two values (most often 0, 1). This is
the response variable and should be of length equal to the number of row of
G1 and G2 arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{n.perm}{integer strictly superior to 0. This is how much permutations
of Y are performed to compute the pvalue. By default, this is fixed to 1000.}
}
\value{
A list with class \code{"GGItest"} containing the following components:

	\item{statistic}{The value of the statistic (minimum p-value, tTS or tProd)}
	\item{p.value}{Tue p-value for the test}
	\item{method}{The method used to perform the Gene-Gene interaction test}
	\item{parameter}{A list of the parameters used to perform the Gene-Gene Interaction test}

}
\description{
\code{GBIGM.test} performs a Gene-Gene Interaction (GGI) analysis by contrasting the information entropy between cases and controls.
}
\details{
First, the conditional entropy and information gain rate is computed for each gene G1 and G2. In a second step, information gain rate for the gene pair (G1,G2) is computed. A p-value is estimated using permutations of Y.
}
\examples{
data(gene.pair)
GBIGM.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2)
}
\references{
	\enumerate{\item J. Li, et al.. A gene-based information gain method for detecting gene-gene interactions in case-control studies. European Journal of Human Genetics, 23 :1566-1572, 2015.
	}
 }

\seealso{
\code{\link{GGI}}
}

