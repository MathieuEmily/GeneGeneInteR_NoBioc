\name{KCCA.test}
\alias{KCCA.test}
\title{KCCA (Kernel Canonical-Correlation Analysis) based GGI analysis.}
\usage{
KCCA.test(Y, G1, G2, kernel=rbfdot(sigma = 0.05), n.boot = 500)
}
\arguments{
\item{Y}{numeric or factor vector with two values (most often 0, 1). This is
the response variable and should be of length equal to the number of row of
G1 and G2 arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{kernel}{A kernel function \code{\link[kernlab]{dots} object. Default is a Gaussian RBF kernel}

\item{n.boot}{integer strictly superior to 0. This is the number of bootstrap
replicates for the variances estimation. By default, this is fixed to 500.}
}
\value{
A list with class \code{"GGItest"} containing the following components:

	\item{statistic}{The value of the statistic (minimum p-value, tTS or tProd)}
	\item{p.value}{Tue p-value for the test}
	\item{method}{The method used to perform the Gene-Gene interaction test}
	\item{parameter}{A list of the parameters used to perform the Gene-Gene Interaction test}
}
\description{
\code{KCCA.test} performs a Gene-Gene Interaction (GGI) analysis based on the
difference of canonical correlations between cases and controls. The "kernel trick" is applied to the canonical correlation to allow the detection of non-linear co-association.
}
\details{
	The test statistic is based on the difference between a Fisher's' transformation of the maximum of the kernelized canonical correlations in cases and controls. To calculate the test statistic for the interaction pvalue, \code{CCA.test} estimates the variance of the Fisher's' transformation of the maximum of the kernelized canonical correlations in cases and controls using a bootstrap method.
	The computation of \code{\link[kernlab]{kcca}.} can be very long.
}
\examples{
\dontrun{
data(gene.pair)
KCCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2)
}
}
\references{
\enumerate{
\item Yuan, Z. et al. (2012): Detection for gene-gene co-association via kernel canonical correlation analysis, BMC Genetics, 13, 83.
\item Larson, N. B. et al. (2013): A kernel regression approach to gene-gene interaction detection for case-control studies, Genetic Epidemiology, 37, 695-703.
}
}
\seealso{
\code{\link{GGI}}, \code{\link{CCA.test}}
}

