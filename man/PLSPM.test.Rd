\name{PLSPM.test}
\alias{PLSPM.test}
\title{PLSPM (Partial Least Squares Path Modelling Analysis) based Gene-Gene Interaction analysis.}
\usage{
PLSPM.test(Y, G1, G2,n.perm=500)
}
\arguments{
\item{Y}{numeric or factor vector with two values (most often 0, 1). This is
the response variable and should be of length equal to the number of row of
G1 and G2 arguments (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{n.perm}{integer strictly superior to 0. This is how much permutations of
Y are performed to compute the pvalue. By default, this is fixed to 500.}
}
\value{
If the function succeed, it returns the interaction pvalue between G1
 and G2. If the pvalue cannot be computed, an error message is displayed.
}
\description{
\code{PLSPM.test} performs a Gene-Gene Interaction (GGI) analysis based on the
modelisation of a statistical relations network. The aim is to quantify the
connections between the latent and the manifest variables.
}
\details{
	The PLSPM-based method, as described in Zhang et al. (2013), aims at comparing the path coefficients \eqn{\beta_D} and \eqn{\beta_C}, where \eqn{\beta_D} is calculated among cases and \eqn{\beta_C} among controls. The co-association between genes \code{G1} and \code{G2} is defined by: \deqn{U=\frac{\beta_D-\beta_C}{\sqrt{Var(\beta_D-\beta_C)}}}.
	The \code{\link[plspm]{plspm}} R package is used to estimate \eqn{U}. The significance pvalue is obtained by using a permutation method on the difference between the path coefficients.
}
\examples{
	data(gene.pair)
PLSPM.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,n.perm=100)
}
\references{
\enumerate{\item X. Zhang et al. (2013) A PLSPM-based test statistic for detecting gene-gene co-association in genome-wide association study with case-control design. PLoS ONE, 8(4) :e62129, April 2013. 
	}
}
\seealso{
\code{\link{GGI}}
}

