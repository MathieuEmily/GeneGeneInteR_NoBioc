\name{CLD.test}
\alias{CLD.test}
\title{CLD (Composite Linkage Desequilibrium) based GGI analysis.}
\usage{
CLD.test(Y, G1, G2, n.perm = 1000)
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
Y are performed to compute the pvalue. By default, this is fixed to 1000.}
}
\value{
If the function succeeds, it returns the interaction pvalue between G1
 and G2. If the pvalue cannot be computed, an error message is displayed. An
 other error message is shown if the test statistic can't be computed.
}
\description{
\code{CLD.test} performs a Gene-Gene Interaction (GGI) analysis based on the
difference between the Composite Linkage Disequilibrium measured in cases and controls respectively.
}
\details{
The test statistic is based on Nagao normalized
Quadratic Distance (NQD). The pvalue is calculated using permutation tests.
}
\examples{
data(gene.pair)
CLD.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2)
}
\references{
\itemize{\item Indika Rajapakse, Michael D. Perlman, Paul J.
 Martin, John A. Hansen, and Charles Kooperberg. Multivariate detection of
 gene-gene interactions. Genetic Epidemiology, 36(6) :622-630, 2012.}
}
\seealso{
\code{\link{GGI}}
}

