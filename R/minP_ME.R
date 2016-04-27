#'@title SNP-based GGI analysis with minP test
#'
#'@description Gene-Gene Interaction analysis is performed based on a SNP-SNP
#'  Interaction (SSI) analysis as a basis. For each gene pairs, a SSI is
#'  performed and one of the following method is used to combine those results
#'  into one GGI p-value: \itemize{ \item Minimum p-value/MinP test \item Gene
#'  Association Test with Extended Simes procedure/ GATES test \item Truncated
#'  Tail Strength/tTS test \item Truncated p-value Product/tProd test}
#'
#'@details A SNP-SNP interaction analysis is first performed between all SNP of
#'  each genes and the results are then summed up into Gene-Gene Interaction
#'  p-values using one of the available methods.
#'
#'  A \eqn{\Sigma} matrix is computed and is used as a multivariate normal
#'  distribution generator. \eqn{\Sigma} is a matrix of size
#'  \eqn{(m1*m2)*(m1*m2)} where \eqn{m1} is the number of SNP in Gene 1 and
#'  \eqn{m2} is the number of SNP in Gene 2. General term \eqn{(i;j)} of
#'  \eqn{\Sigma} is \eqn{m1:m2[i] * m1:m2[j]} where m1:m2 is the ordered vector
#'  of linkage desequilibrium (LD) R coefficient bewteen all SNP of the two
#'  genes. Corresponding multivariate normal distribution, \eqn{\Phi}, is
#'  integrated.
#'
#'  \itemize{ \item Minimum p-value Test (MinP test) - MinP test is based on the
#'  minimum SNP-SNP interaction p-value observed. MinP p-value is calculated by
#'  integrating \eqn{\Phi}. \item Gene Association Test using Extended Simes
#'  procedure (GATES test) - As with MinP, only strongest signal is considered.
#'  Strongest signal is not defined as minimum p-value but as the minimum value
#'  of \deqn{me*p[i]/m[i]} where me is the effective number of independant
#'  tests, p[i] is the i-th p-values and m[i] is the effective number of
#'  independant test among the top i p-values. P-values are ascending ordered.
#'  GATES p-value is thus the minimum term \eqn{me*p[i]/m[i]}. Multiple methods
#'  exists to estimate \eqn{me} and \eqn{m[i]} terms: \itemize{\item
#'  Cheverud-Nyholt method \item Keff method \item Li & Ji method} See
#'  references to get more details about those methods. \item Truncated Tail
#'  Strenght Test (tTS test) - tTS test does not consider only the strongest
#'  signal but all signals that are inferior to a threshold \eqn{\tau}. For these
#'  p-values, the weighted sum of \deqn{1-p[i]*m1*m2/i} is computed and
#'  represents the test statistic. The empirical p-value is calculated using
#'  \code{\link[mvtnorm]{rmvnorm}} with the sigma matrix. \item Truncated
#'  p-value Product Test (TProd test) - This is the same as tTS but the test
#'  statistic is a product of the p-values, not a sum.}
#'
#'@param Y numeric or factor vector with two values (most often 0, 1). This is
#'  the response variable and should be of length equal to the number of row of
#'  G1 and G2 arguments (number of individuals).
#'@param G1,G2 \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix}
#'  objects. Must have a number of row equal to Y argument's length.
#'@param me.est (optional, \emph{specific to \code{gates.test}}) character
#'  string for GATES method. Me parameter estimation method. Must be one of the
#'  following: "ChevNy", "Keff", "LiJi" (See details)
#'@param alpha \emph{(optional, specific to \code{gates.test})} numeric value in [0, 1].
#'  Threshold for GATES method when estimating Me with Keff method.
#'@param tau \emph{(optional, specific to \code{tTS.test} and
#'  \code{tProd.test})} numeric in [0, 1]. See details section for its use.
#'@param n.sim \emph{(optional, specific to \code{tTS.test} and
#'  \code{tProd.test})} numeric positive value. Number of observations for the
#'  multivariate normal distribution simulation to compute p-values for tTS and
#'  tProd methods.
#'
#'@return GGI analysis p-value is returned.
#'
#'@references Ma L., Clark AG., Keinan A. (2013) Gene-Based Testing Of
#'  Interactions in Association Studies of Quantitative Traits. PLoS Genet
#'  9(2):e1003321.doi: 10.1371/journal.pgen.1003321
#'
#'  Moskvina V., Schmidt KM. (Sept. 2008) On multiple-testing correction in
#'  genome-wide association studies. Genet Epidemiol., 32(6):567-73. doi:
#'  10.1002/gepi.20331.
#'
#'  J. Li, L. Ji, (Aug. 2005) Adjusting multiple testing in multilocus analyses
#'  using the eigenvalues of a correlation matrix. Heredity 95, 221-227
#'
#'@seealso \code{\link{GGI}}, \code{\link{GGI.plot}}
#'
#'@export
#'
#'@examples
#' ## Dataset is included in the package
#' ped <- system.file("extdata/example.ped", package="GGItest")
#' info <- system.file("extdata/example.info", package="GGItest")
#' posi <- system.file("extdata/example.txt", package="GGItest")
#' dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\t")
#'
#' G1 <- select.snps(dta$snpX, dta$genes.info, "TXNDC5")$snpX
#' G2 <- select.snps(dta$snpX, dta$genes.info, "DNAH9")$snpX
#'
#' resp <- system.file("extdata/response.txt", package="GGItest")
#' Y  <- read.csv(resp, header=F)
#' Y <- Y[,1]
#'
#' ## All SSI methods are very similar in their use
#' minP.test(Y, G1, G2)
#' gates.test(Y, G1, G2, me.est = "LiJi", alpha = 0.2)
#' tTS.test(Y, G1, G2, tau = 0.5, n.sim = 500)
#' tProd.test(Y, G1, G2, tau = 0.5, n.sim = 500)
minP.test <- function(Y, G1, G2){

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  if (nlevels(as.factor(Y)) != 2) {
    stop("response variable should be binary. (2 modes).")
  } else if (class(G1) != "SnpMatrix" | class(G2) != "SnpMatrix") {
    stop("G1 and G2 arguments should be SnpMatrix objects.")
  } else if (nrow(G1) != nrow(G2)) {
    stop("G1 and G2 should have same rows count.")
  } else if (length(Y) != nrow(G1) | length(Y) != nrow(G2)) {
    stop("Response variable should be conformant with genes matrices rows number.")
  } else if (sum(is.na(G1))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(G2))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(Y)) != 0) {
    stop("The response variable vector must be complete. No NAs are allowed.")
  }

  Y <- as.numeric(Y)
  if(min(Y)!=0){Y<-Y-min(Y)}
  Y <- as.factor(Y)

  # minP test is set up
  # SNP interactions are computed
  SSI.res <- SSI.test(Y, G1, G2)

  ## Computation of the correlation matrix
	MatCor1 <- 	snpStats::ld(G1, G1, stats="R")
	MatCor2 <- 	snpStats::ld(G2, G2, stats="R")
	n1 <- ncol(MatCor1)
	n2 <- ncol(MatCor2)
	n.pairs <- n1*n2
	sigma.matrix <- matrix(NA,ncol=n.pairs,nrow=n.pairs)
	for (i in 1:(n.pairs-1)){
		i1 <- floor((i-1)/n2)+1
		j1 <- i-(i1-1)*n2
		for (j in (i+1):n.pairs){
			i2 <- floor((j-1)/n2)+1
			j2 <- j-(i2-1)*n2
			sigma.matrix[i,j] <- sigma.matrix[j,i] <- MatCor1[i1,i2]*MatCor2[j1,j2]
		}
	}
	diag(sigma.matrix) <- 1

  # minP test is computed past this point
  SSI.min <- min(SSI.res, na.rm=T)

  # Case in which 1 - SSI.min/2 == 1 (qnorm(1) = -Inf)
  if (SSI.min == 0) {
    GG.Pmin <- 0
    # Case in which 1 - SSI.min/2 == 1 (qnorm(1) = Inf)
  } else if (SSI.min > 0 & SSI.min < 0.5) {
    lower <- rep(qnorm(SSI.min / 2), ncol(sigma.matrix))
    upper <- rep(qnorm(1 - SSI.min / 2), ncol(sigma.matrix))

    GG.Pmin <- 1 - mvtnorm::pmvnorm(lower=lower,upper=upper,sigma=sigma.matrix,maxpts=2.5E6,abseps = 0.01
    #abseps=1E-13
    )
  }
  if (attributes(GG.Pmin)$msg=="Completion with error > abseps"){
		warning(paste("p-values are approximated with error=",attributes(GG.Pmin)$error))
	}
	pval <- as.numeric(GG.Pmin)
	stat <- SSI.min
	names(stat)="minP"
	res <- list(statistic=stat,p.value=pval,method="minP")
	class(res) <- "GGItest"
  return(res)

}
