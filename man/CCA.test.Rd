\name{CCA.test}
\alias{CCA.test}
\title{CCA (Canonical-Correlation Analysis) based GGI analysis.}
\usage{
CCA.test(Y, G1, G2, n.boot = 500)
}
\arguments{
\item{Y}{numeric or factor vector with two values (most often 0, 1). This is
the response variable and should be of length equal to the number of row of
G1 and G2 arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{n.boot}{integer strictly superior to 0. This is the number of bootstrap
replicates for estimating variances. By default, this is fixed to 500.}
}
\value{
If the function succeed, it returns the interaction pvalue between G1
 and G2. An error message is displayed if the variances used for the test
 statistic calculus could not be estimated. An error message is shown if the
 test statistic could not be computed.
}
\description{
\code{CCA.test} performs a Gene-Gene Interaction (GGI) analysis based on the
difference of canonical correlation between cases and controls.
}
\details{
	The test statistic is based on the difference between a Fisher's' transformation of the maximum of the canonical correlations in cases and controls. To calculate the test statistic for the interaction pvalue, \code{CCA.test} estimates the variance of the Fisher's' transformation of the maximum of the canonical correlations in cases and controls using a bootstrap method.
}
\examples{
	data(gene.pair)
	CCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2)
}
\references{
\enumerate{\item Qianqian Peng, Jinghua Zhao, and Fuzhong Xue. A
 gene-based method for detecting gene-gene co-association in a case-control
 study. European Journal of Human Genetics, 18(5) :582-587, May 2010.}
}
\seealso{
\code{\link{GGI}}, \code{\link{KCCA.test}}
}

