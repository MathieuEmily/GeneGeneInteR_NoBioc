#'PLSPM (Partial Least Squares Path Modelling Analysis) based GGI analysis.
#'
#'\code{PLSPM.test} performs a Gene-Gene Interaction (GGI) analysis based on the
#'modelisation of a statistical relations network. The aim is to quantify the
#'connections between the latent and the manifest variables.
#'
#'To calculate the test statistic for the interaction pvalue, \code{PLSPM.test}
#'uses \code{\link[plspm]{plspm}} for path coefficients and variances
#'estimation. The pvalue is obtained by using a permutation method on the
#'difference between the path coefficients.
#'
#'@param Y numeric or factor vector with two values (most often 0, 1). This is
#'  the response variable and should be of length equal to the number of row of
#'  G1 and G2 arguments (number of individuals).
#'@param G1 \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
#'  Must have a number of row equal to Y argument's length.
#'@param G2 \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
#'  Must have a number of row equal to Y argument's length.
#'
#'@return If the function succeed, it returns the interaction pvalue between G1
#'  and G2. If the pvalue cannot be computed, an error message is displayed.
#'
#'@references \enumerate{\item Xiaoshuai Zhang, Xiaowei Yang, Zhongshang Yuan,
#'  Yanxun Liu, Fangyu Li, Bin Peng, Dianwen Zhu, Jinghua Zhao, and Fuzhong Xue.
#'  A plspm-based test statistic for detecting gene-gene co-association in
#'  genome-wide association study with case-control design. PLoS ONE, 8(4)
#'  :e62129, April 2013. \item H. Wold. Estimation of principal components and
#'  related models by iterative least squares. Academic Press, pages 391-420,
#'  1966. \item Doerge RW and Churchill GA. Permutation tests for multiple loci
#'  affecting a quantitative character. Genetics, (142) :285, 1996. \item M.
#'  Hallin and C. Ley. Permutation tests. Wiley StatsRef : Statistics Reference
#'  Online, 2014.}
#'
#'@seealso \code{\link{GGI}}
#'
#'@export
#'
#' @examples
#' ## Dataset is included in the package
#' ped <- system.file("extdata/example.ped", package="GGItest")
#' info <- system.file("extdata/example.info", package="GGItest")
#' posi <- system.file("extdata/example.txt", package="GGItest")
#' dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\t")
#'
#' dta <- imputeSnpMatrix(dta$snpX, genes.info = dta$genes.info)
#'
#' G1 <- select.snps(dta$snpX, dta$genes.info, "bub3")$snpX
#' G2 <- select.snps(dta$snpX, dta$genes.info, "CA1")$snpX
#'
#' resp <- system.file("extdata/response.txt", package="GGItest")
#' Y  <- read.csv(resp, header=FALSE)
#'
#' ## By default, the number of bootstrap replicates is 500
#' PLSPM.test(Y, G1, G2)

PLSPM.test <- function(Y, G1, G2, n.perm=500){

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  if(nlevels(as.factor(Y))!=2){
    stop("Y must be a factor with 2 levels, most likely 0 and 1.")
  } else if(class(G1)!="SnpMatrix"){
    stop("G1 must be a SnpMatrix object.")
  } else if(class(G2)!="SnpMatrix"){
    stop("G2 must be a SnpMatrix object")
  } else if(nrow(G1)!=nrow(G2)){
    stop("Both G1 and G2 must contain the same number of individuals.")
  } else if(length(Y)!=nrow(G1)){
    stop("Y and both SnpMatrix objects must contain the same number of individuals.")
  } else if (sum(is.na(G1))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(G2))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (sum(is.na(Y))!=0) {
    stop("The response variable must be complete. No NAs are allowed.")
  }

  X1 <- as(G1, "numeric")
  X2 <- as(G2, "numeric")

  Y <- as.numeric(Y)
  if(min(Y)!=0){Y<-Y-min(Y)}
  Y <- as.factor(Y)

  X <- cbind(X1,X2)

  Gene1 <- c(0,0)
  Gene2 <- c(1,0)

  w1 <- which(Y==1)
  w0 <- which(Y==0)

  XCases <- X[w1,]
  XControls <- X[w0,]

  my.path <- rbind(Gene1,Gene2)
  my.blocks <- list(1:ncol(X1),(ncol(X1)+1):(ncol(X1)+ncol(X2)))
  my.modes = c("A", "A")

  mod1<-NULL;
  mod0<-NULL;

  try(mod1 <- plspm::plspm(XCases,my.path,my.blocks, modes = my.modes), silent=T)
  if(is.null(mod1)){warning("P-value could not be computed. NA returned");return(NA)}

  try(mod0 <- plspm::plspm(XControls,my.path,my.blocks, modes = my.modes),silent=T)
  if(is.null(mod0)){warning("P-value could not be computed. NA returned");return(NA)}

  beta1 <- mod1$inner_model[[1]][2,1]
  vbeta1 <- mod1$inner_model[[1]][2,2]^2
  beta0 <- mod0$inner_model[[1]][2,1]
  vbeta0 <- mod0$inner_model[[1]][2,2]^2

  U <- (beta0-beta1)/sqrt(vbeta0+vbeta1)
  
  U.perm <- rep(NA,times=n.perm)
  for (i in 1:n.perm){
    restart<-T
    while(restart){
	  	Y.perm <- sample(Y)
  		w1 <- which(Y.perm==1)
	  	w0 <- which(Y.perm==0)
  		XCases <- X[w1,]
	  	XControls <- X[w0,]
  		mod1<-NULL;
		mod0<-NULL;
		try(mod1 <- plspm::plspm(XCases,my.path,my.blocks, modes = my.modes), silent=T)
#		if(is.null(mod1)){warning("P-value could not be computed. NA returned");return(NA)}
		try(mod0 <- plspm::plspm(XControls,my.path,my.blocks, modes = my.modes),silent=T)
#		if(is.null(mod0)){warning("P-value could not be computed. NA returned");return(NA)}
		if (!is.null(mod1) & !is.null(mod0)){restart <- F}
	}
	beta1 <- mod1$inner_model[[1]][2,1]
	vbeta1 <- mod1$inner_model[[1]][2,2]^2
	beta0 <- mod0$inner_model[[1]][2,1]
	vbeta0 <- mod0$inner_model[[1]][2,2]^2
	
	U.perm[i] <- (beta0-beta1)/sqrt(vbeta0+vbeta1)
  }
  
  #pval <- 2*(1-pnorm(abs(U)))
	  pval <- mean(abs(U.perm) > abs(U))
	  stat <- U
	names(stat)="U"
	res <- list(statistic=stat,p.value=pval,method="Partial Least Squares Path Modeling")
	class(res) <- "GGItest"
  return(res)
#  return(pval)
}
