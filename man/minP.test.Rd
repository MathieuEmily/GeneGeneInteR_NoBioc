\name{gates.test}
\alias{gates.test}
\title{SNP-based GGI analysis with minP test}
\usage{
	
	
gates.test(Y, G1, G2, alpha = 0.05, me.est = c("ChevNy", "Keff", "LiJi","Galwey"))

minP.test(Y, G1, G2)

tProd.test(Y, G1, G2, tau = 0.05, n.sim = 1000)

tTS.test(Y, G1, G2, tau = 0.05, n.sim = 1000)
}
\arguments{
\item{Y}{numeric or factor vector with two values (most often 0, 1). This is
the response variable and should be of length equal to the number of row of
G1 and G2 arguments (number of individuals).}

\item{G1, G2}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix}
objects. Must have a number of row equal to Y argument's length.}

\item{alpha}{\emph{(optional, specific to \code{gates.test})} numeric value in [0, 1].
Threshold for GATES method when estimating Me with Keff method.}

\item{me.est}{(optional, \emph{specific to \code{gates.test}}) character
string for GATES method. Me parameter estimation method. Must be one of the
following: "ChevNy", "Keff", "LiJi", "Galwey" (See details)}

\item{tau}{\emph{(optional, specific to \code{tTS.test} and
\code{tProd.test})} numeric in [0, 1]. See details section for its use.}

\item{n.sim}{\emph{(optional, specific to \code{tTS.test} and
\code{tProd.test})} numeric positive value. Number of observations for the
multivariate normal distribution simulation to compute p-values for tTS and
tProd methods.}
}
\value{
GGI analysis p-value is returned.
}
\description{
Gene-Gene Interaction analysis is performed based on a SNP-SNP
 Interaction (SSI) analysis as a basis. For each gene pairs, a SSI is
 performed and one of the following method is used to combine those results
 into one GGI p-value: \itemize{ \item Minimum p-value/MinP test \item Gene
 Association Test with Extended Simes procedure/ GATES test \item Truncated
 Tail Strength/tTS test \item Truncated p-value Product/tProd test}
}
\details{
A SNP-SNP interaction analysis is first performed between all SNP of
 each genes and the results are then summed up into Gene-Gene Interaction
 p-values using one of the available methods.

 A \eqn{\Sigma} matrix is computed and is used as a multivariate normal
 distribution generator. \eqn{\Sigma} is a matrix of size
 \eqn{(m1*m2)*(m1*m2)} where \eqn{m1} is the number of SNP in Gene 1 and
 \eqn{m2} is the number of SNP in Gene 2. General term \eqn{(i;j)} of
 \eqn{\Sigma} is \eqn{m1:m2[i] * m1:m2[j]} where m1:m2 is the ordered vector
 of linkage desequilibrium (LD) R coefficient bewteen all SNP of the two
 genes. Corresponding multivariate normal distribution, \eqn{\Phi}, is
 integrated.

 \itemize{ \item Minimum p-value Test (MinP test) - MinP test is based on the
 minimum SNP-SNP interaction p-value observed. MinP p-value is calculated by
 integrating \eqn{\Phi}. \item Gene Association Test using Extended Simes
 procedure (GATES test) - As with MinP, only strongest signal is considered.
 Strongest signal is not defined as minimum p-value but as the minimum value
 of \deqn{me*p[i]/m[i]} where me is the effective number of independant
 tests, p[i] is the i-th p-values and m[i] is the effective number of
 independant test among the top i p-values. P-values are ascending ordered.
 GATES p-value is thus the minimum term \eqn{me*p[i]/m[i]}. Multiple methods
 exists to estimate \eqn{me} and \eqn{m[i]} terms: \itemize{\item
 Cheverud-Nyholt method \item Keff method \item Li & Ji method} See
 references to get more details about those methods. \item Truncated Tail
 Strenght Test (tTS test) - tTS test does not consider only the strongest
 signal but all signals that are inferior to a threshold \eqn{\tau}. For these
 p-values, the weighted sum of \deqn{1-p[i]*m1*m2/i} is computed and
 represents the test statistic. The empirical p-value is calculated using
 \code{\link[mvtnorm]{rmvnorm}} with the sigma matrix. \item Truncated
 p-value Product Test (TProd test) - This is the same as tTS but the test
 statistic is a product of the p-values, not a sum.}
}
\section{Functions}{
\itemize{
\item \code{gates.test}: SNP-based GGI analysis with GATES test

\item \code{tProd.test}: SNP-based GGI analysis with tProd test

\item \code{tTS.test}: SNP-based GGI analysis with tTS test
}}
\examples{
data(gene.pair)

## All SSI methods are very similar in their use
minP.test(gene.pair$Y, gene.pair$G1, gene.pair$G2)
gates.test(gene.pair$Y, gene.pair$G1, gene.pair$G2, me.est = "Keff",alpha=0.05)
tTS.test(gene.pair$Y, gene.pair$G1, gene.pair$G2, tau = 0.5, n.sim = 500)
tProd.test(gene.pair$Y, gene.pair$G1, gene.pair$G2, tau = 0.5, n.sim = 500)
}
\references{
\enumerate{
	\item M. Emily AGGrEGATOr: A Gene-based GEne-Gene interActTiOn test for case-control association studies, Statistical Application in Genetics and Molecular Biology, 15(2): 151-171, 2016.
	\item L. Ma, A.G. Clark and A. Keinan Gene-Based Testing Of Interactions in Association Studies of Quantitative Traits. PLoS Genetics 9(2):e1003321, 2013.
	\item V. Moskvina and K.M. Schmidt On multiple-testing correction in genome-wide association studies. Genetic Epidemiology, 32(6): 567-73, 2008.
  \item J. Li and L. Ji. Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix. Heredity 95: 221-227, 2005.
  \item N.W. Galwey. A new measure of the effective number of tests, a practical tool for comparing families of non-independent significance tests. Genetic Epidemiology 33(7): 559-568, 2009.
\item J.M. Cheverud JM. A simple correction for multiple comparisons in interval mapping genome scans. Heredity. 87(1):52–58, 2001.
\item D.R. Nyholt. A Simple Correction for Multiple Testing for Single-Nucleotide Polymorphisms in Linkage Disequilibrium with Each Other. American journal of human genetics. 74(4): 765–769, 2004.
 }
}
\seealso{
\code{\link{GGI}}
}

