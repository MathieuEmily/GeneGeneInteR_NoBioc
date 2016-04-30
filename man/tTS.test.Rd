\name{tTS.test}
\alias{tTS.test}
\title{Gene-based Gene-Gene Interaction analysis by combining SNP-SNP interaction tests}
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
Threshold for GATES method when estimating the number of effective tests Me with Keff method.}

\item{me.est}{(optional, \emph{specific to \code{gates.test}}) character
string for GATES method. \code{me.est} corresponds to the method for estimating the number of effective tests. Must be one of the
following: "ChevNy", "Keff", "LiJi", "Galwey" (See details)}

\item{tau}{\emph{(optional, specific to \code{tTS.test} and
\code{tProd.test})} numeric in [0, 1]. See details section for its use.}

\item{n.sim}{\emph{(optional, specific to \code{tTS.test} and
\code{tProd.test})} numeric positive value. \code{n.sim} is the number of multivariate normal distribution simulations used to compute the significant level for tTS and
tProd methods.}
}
\value{
GGI analysis p-value is returned.
}
\description{
\code{gates.test}, \code{minP.test}, \code{tTS.test} and \code{tProd.test} aim at performing gene-gene interaction analysis based on SNP-SNP interaction tests. The following methods are used to combine SNP-SNP interaction tests
 into a single Gene-Gene Interaction p-value: 
 \itemize{ 
 	\item Minimum p-value in \code{minP.test} function 
 	\item Gene Association Test with Extended Simes procedure in \code{gates.test}
 	\item Truncated Tail Strength in \code{tTS.test} function
 	\item Truncated p-value Product in \code{tProd.test} function
}
}
\details{
	In a first step, all methods start by applying a logistic regression model to test all pairs of SNPs between the two genes \code{G1} and \code{G2}. If \code{G1} has \eqn{m1} SNPs and \code{G2} \eqn{m2} SNPs, a total of \eqn{m1 * m2} SNP-SNP tests are performed. In a second step, the \eqn{m1 * m2} SNP-SNP tests are combined according to their covariance matrix \eqn{\Sigma}. \eqn{\Sigma} is computed as described in the method developped in Emily (2016). The covariance \eqn{\Sigma} is used in each method as follows:

 \itemize{ 
 	\item minP test - minP test considered the significant of the observed minimum p-value. Significance is computed by integrating the multivariate normal distribution with covariance \eqn{\Sigma} as proposed in Conneelly and Boehnke (2008).
 	\item GATES test - The p-value for GATES is the minimum p-value obtained after a multiple testing correction of the SNP-SNP interaction p-value. Correction for multiple testing is defined as \deqn{me * p[i]/me[i]} where me is the effective number of independant tests, p[i] is the i-th top p-values and me[i] is the effective number of
 independant test among the top i p-values. Many methods
 exist to estimate \eqn{me} and \eqn{me[i]} terms: 
 \itemize{
 	\item Cheverud-Nyholt method (Cheverud, 2001 and Nyholt, 2004)
 	\item Keff method (Moskovina and Schmidt, 2008)
 	\item Li & Ji method (Li and Ji, 2005)
 	\item Galwey method (Galwey, 2009)
 } Details of each method can be found in the references. 
 	\item tTS test - tTS test does not consider only the strongest signal but all signals that are inferior to a given threshold \eqn{\tau}. For these
 p-values, the weighted sum of \deqn{tTS=\sum (1-p[i]*(m1*m2+1)/i)} is computed and
 represents the test statistic. The p-value is calculated using an empirical distribution of \eqn{tTS} obtained by simulating multivariate normal statistics with a covariance \eqn{\Sigma} as proposed by Jiang et al. (2011). 
 \item TProd test - The procedure is similar to \eqn{tTS} with \eqn{tProd} defined as \deqn{tProd=\prod p[i]}
. See Zaykin et al. (2002) for details.}
}
\examples{
data(gene.pair)

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
\item J.M. Cheverud. A simple correction for multiple comparisons in interval mapping genome scans. Heredity. 87(1):52–58, 2001.
\item D.R. Nyholt. A Simple Correction for Multiple Testing for Single-Nucleotide Polymorphisms in Linkage Disequilibrium with Each Other. American journal of human genetics. 74(4): 765–769, 2004.
\item K.N. Conneely and M. Boehnke. So many correlated tests, so little time! rapid adjustment of p values for multiple correlated tests. The American Journal of Human Genetics, 81: 1158–1168, 2008
\item B. Jiang, X. Zhang, Y. Zuo and G. Kang. A powerful truncated tail strength method for testing multiple null hypotheses in one dataset. Journal of Theoretical Biology 277: 67–73, 2011.
\item D.V. Zaykin, L.A. Zhivotovsky, P.H. Westfall and B.S. Weir. Truncated product method for combining P-values. Genetic epidemiology 22: 170–185, 2002.
 }
}
\seealso{
\code{\link{GGI}}
}

