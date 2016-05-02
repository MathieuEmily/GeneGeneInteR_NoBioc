\name{PCA.GenFreq}
\alias{PCA.GenFreq}
\title{PCA-based Gene-Gene Interaction analysis.}
\usage{
PCA.Std(Y, G1, G2, threshold = 0.8)

PCA.GenFreq(Y, G1, G2, threshold = 0.8)
}
\arguments{
\item{Y}{numeric or factor vector with two values (most often 0, 1). This is
the response variable and should be of length equal to the number of row of
\code{G1} and \code{G2} (number of individuals).}

\item{G1}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
Must have a number of row equal to Y argument's length.}

\item{threshold}{(optional) numeric comprised in \eqn{[0, 1]} interval. This is how much inertia should be kept for each PCA and has a direct impact on how many PCA dimensions are kept.}
}
\value{
If the function succeed, the interaction p-value is returned.

}
\description{
\code{PCA.Std} \code{PCA.GenFreq} perform a Gene-Gene Interaction (GGI) analysis by testing the interaction between the principal components of the two genes.
 using Principal Components AnalysisPCA in which dataset is standardized using variables' standard
 deviation. With \code{PCA.GenFreq}, dataset is standardized using the
 standard deviation under Hardy-Weinberg equilibrium.
}
\details{
In a first step, \code{PCA.Std} and \code{PCA.GenFreq} performs a Principal Components Analysis on both \code{G1} and \code{G2} genes that are seen as matrices of allele counts. With \code{PCA.Std}, dataset is standardized using variables' standard deviation,  while dataset is standardized using standard deviation under Hardy-Weinberg equilibrium for  \code{PCA.GenFreq}. 

 Principal components are retrieved to describe each dataset with user-set
 inertia percentage (parameter \code{threshold}).

 Principal components are then used in a logit model fitting process in which
 interaction is tested. If model fitting fails then a PC is removed from the
 pool and a new attempt to fitting is made. The process loops until either
 the model is fitted or there are only two PC for each genes remaining.
 Principal components are removed according to lowest inertia percentage.
}
\examples{
data(gene.pair)
PCA.Std(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,threshold=0.7)
PCA.GenFreq(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,threshold=0.7)
}
\references{
J. Li et al. (2009) Identification of gene-gene interaction using principal components. BMC
 Proceedings, 3 (Suppl. 7): S78
}
\seealso{
\code{\link{GGI}}
}

