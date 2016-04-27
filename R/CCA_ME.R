#'CCA (Canonical-Correlation Analysis) based GGI analysis.
#'
#'\code{CCA.test} performs a Gene-Gene Interaction (GGI) analysis based on the
#'difference of correlation between case and control.
#'
#'To calculate the test statistic for the interaction pvalue, \code{CCA.test}
#'estimate the variances of each group case-control. A bootstrap method is used
#'to compute these variances. With the difference between the maximum of the
#'canonical correlations for each group tranformed using an analogy to Fisher's
#'transformation, the test statistic is computed and represents the
#'co-association between genes.
#'
#'@param Y numeric or factor vector with two values (most often 0, 1). This is
#'  the response variable and should be of length equal to the number of row of
#'  G1 and G2 arguments (number of individuals).
#'@param G1 \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
#'  Must have a number of row equal to Y argument's length.
#'@param G2 \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
#'  Must have a number of row equal to Y argument's length.
#'@param n.boot integer strictly superior to 0. This is the number of bootstrap
#'  replicates for the variances estimation. By default, this is fixed to 500.
#'
#'@return If the function succeed, it returns the interaction pvalue between G1
#'  and G2. An error message is displayed if the variances used for the test
#'  statistic calculus could not be estimated. An error message is shown if the
#'  test statistic could not be computed.
#'
#'@references \enumerate{\item Qianqian Peng, Jinghua Zhao, and Fuzhong Xue. A
#'  gene-based method for detecting gene-gene co-association in a case-control
#'  study. European Journal of Human Genetics, 18(5) :582-587, May 2010. \item
#'  R. A. Fisher. On the probable error of a coefficient of correlation deduced
#'  from a small sample. Metron, 1 :3-32, 1921.}
#'
#'@seealso \code{\link{GGI}}, \code{\link{KCCA.test}}
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
#' G1 <- select.snps(dta$snpX, dta$genes.info, "TXNDC5")$snpX
#' G2 <- select.snps(dta$snpX, dta$genes.info, "DNAH9")$snpX
#'
#' resp <- system.file("extdata/response.txt", package="GGItest")
#' Y  <- read.csv(resp, header=FALSE)
#'
#' ## By default, the number of bootstrap replicates is 500
#' CCA.test(Y, G1, G2)

CCA.test <- function(Y, G1, G2, n.boot = 500){

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  if(class(n.boot)!="numeric"){
    stop("n.boot must be numeric.")
  } else if(n.boot<1){
    stop("n.boot must be strictly superior to 0.")
  } else if(nlevels(as.factor(Y))!=2){
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

  stat<-"empty"
  try(stat <- get.U(Y=Y,X1=X1,X2=X2,n.boot=n.boot))

  if(is.na(stat)){
    pval<-1

  } else {
    if(stat=="empty"){stop("P-value could not be computed, test statistic is missing")}
    pval <- 2*(1-pnorm(abs(stat)))

  }
	names(stat)="CCU"
	res <- list(statistic=stat,p.value=pval,method="Canonical Correlation Analysis")
	class(res) <- "GGItest"
  return(res)
}


get.U <- function(Y,X1,X2,n.boot=500){
  w0 <- which(Y==0)
  w1 <- which(Y==1)
  X1.0 <- X1[w0,]
  X2.0 <- X2[w0,]
  X1.1 <- X1[w1,]
  X2.1 <- X2[w1,]
  z0 <- get.z(X1.0,X2.0)
  z1 <- get.z(X1.1,X2.1)
  if(is.na(z0)&&is.na(z1)){
    warning("Canonical correlations between the first gene cases and the second gene cases equals to 1")
    warning("Canonical correlations between the first gene controls and the second gene controls equals to 1")
    return(0)
  }else{
    if(is.na(z0)){
      warning("Canonical correlations between the first gene controls and the second gene controls equals to 1")
      return(NA)
    }else{
      if(is.na(z1)){
        warning("Canonical correlations between the first gene cases and the second gene cases equals to 1")
        return(NA)
      }else{
        vz0<-NA;vz1<-NA;
        try(vz0 <- estim.var.z(X1.0,X2.0,n.boot=n.boot))
        try(vz1 <- estim.var.z(X1.1,X2.1,n.boot=n.boot))
        if(is.na(vz0)||is.na(vz1)){stop("The test statistic could not be computed, the variance estimator is missing")}
        else{return((z0-z1)/sqrt(vz0+vz1))}
      }
    }
  }
}


get.z <- function(X1,X2){
  tmp <- cancor(X1,X2)
  r <- tmp$cor[1]
  if(r>1-10^-16){
    z=NA;
  }else{
    z <- (1/2)*(log(1+r)-log(1-r))
  }
  return(z)
}

estim.var.z <- function(X1,X2,n.boot=500){
  z.vec <- rep(NA,times=n.boot)
  for (i in 1:n.boot){
    restart<-T
    while(restart){
      ind.boot <- sample(1:nrow(X1),nrow(X1),replace=TRUE)
      X1.boot <- as.matrix(X1[ind.boot,])
      X2.boot <- as.matrix(X2[ind.boot,])
      z.vec[i] <- get.z(X1.boot,X2.boot)
      if(!is.na(z.vec[i])){restart<-F}
    }
  }
  return(var(z.vec))
}
