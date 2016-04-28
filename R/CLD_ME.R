#'CLD (Composite Linkage Desequilibrium) based GGI analysis.
#'
#'\code{CLD.test} performs a Gene-Gene Interaction (GGI) analysis based on the
#'difference between the covariance matrices of case and control.
#'
#'The test statistic is computed using a method based on Nagao normalized
#'Quadratic Distance (NQD). The pvalue is calculated by bootstrap method.
#'
#'@param Y numeric or factor vector with two values (most often 0, 1). This is
#'  the response variable and should be of length equal to the number of row of
#'  G1 and G2 arguments (number of individuals).
#'@param G1 \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
#'  Must have a number of row equal to Y argument's length.
#'@param G2 \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
#'  Must have a number of row equal to Y argument's length.
#'@param n.perm integer strictly superior to 0. This is how much permutations of
#'  Y are performed to compute the pvalue. By default, this is fixed to 1000.
#'
#'@return If the function succeed, it returns the interaction pvalue between G1
#'  and G2. If the pvalue cannot be computed, an error message is displayed. An
#'  other error message is shown if the test statistic can't be computed.
#'
#'@references \enumerate{\item Indika Rajapakse, Michael D. Perlman, Paul J.
#'  Martin, John A. Hansen, and Charles Kooperberg. Multivariate detection of
#'  gene-gene interactions. Genetic Epidemiology, 36(6) :622-630, 2012. \item
#'  Hisao Nagao. On some test criteria for covariance matrix. Ann. Statist.,
#'  1(4) :700-709, July 1973.}
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
#' G1 <- select.snps(dta$snpX, dta$genes.info, "TXNDC5")$snpX
#' G2 <- select.snps(dta$snpX, dta$genes.info, "DNAH9")$snpX
#'
#' resp <- system.file("extdata/response.txt", package="GGItest")
#' Y  <- read.csv(resp, header=FALSE)
#'
#' ## By default, the number of bootstrap replicates is 500
#' CLD.test(Y, G1, G2)

## toctoc

CLD.test <- function(Y, G1, G2, n.perm = 1000){

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  if(class(n.perm)!="numeric"){
    stop("n.boot must be numeric.")
  } else if(n.perm<1){
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

  stat <-NA
  try(stat <- CLD(Y,X1,X2,0))
  if(is.na(stat)){stop("P-value can't be computed")}

  stat.perm <- rep(NA,times=n.perm)
  for (i in 1:n.perm){
    while(is.na(stat.perm[i])){
      Y.perm <- sample(Y,length(Y))
      stat.perm[i] <- CLD(Y.perm,X1,X2,1)
    }
  }
  
  pval <- mean(stat.perm > stat)
  	names(stat)="CLD"
	res <- list(statistic=stat,p.value=pval,method="Composite Linkage Disequilibrium")
	class(res) <- "GGItest"
  return(res)
}

CLD <- function(Y,X1,X2,bool){
  p1 <- ncol(X1)
  p2 <- ncol(X2)

  wY1 <- which(Y==1)
  n <- length(wY1)
  wY0 <- which(Y==0)
  m <- length(wY0)

  X1.1 <- X1[wY1,]
  X1.0 <- X1[wY0,]
  X2.1 <- X2[wY1,]
  X2.0 <- X2[wY0,]

  S <- cov(cbind(X1.1/2,X2.1/2))
  S11 <- S[1:p1,1:p1]
  S22 <- S[(p1+1):(p1+p2),(p1+1):(p1+p2)]
  S12 <- S[(p1+1):(p1+p2),1:p1]
  S21 <- S[1:p1,(p1+1):(p1+p2)]
  my.T <- cov(cbind(X1.0/2,X2.0/2))
  T11 <- my.T[1:p1,1:p1]
  T22 <- my.T[(p1+1):(p1+p2),(p1+1):(p1+p2)]
  T12 <- my.T[(p1+1):(p1+p2),1:p1]
  T21 <- my.T[1:p1,(p1+1):(p1+p2)]

  W <- (m*S+n*my.T)/(m+n)
  W11 <- W[1:p1,1:p1]
  W22 <- W[(p1+1):(p1+p2),(p1+1):(p1+p2)]

  STilde <- W
  STilde[(p1+1):(p1+p2),1:p1]  <- S12
  STilde[1:p1,(p1+1):(p1+p2)]  <- S21

  TTilde <- W
  TTilde[(p1+1):(p1+p2),1:p1]  <- T12
  TTilde[1:p1,(p1+1):(p1+p2)]  <- T21

  delta2<-NA
  if(!bool){
    try(delta2 <- sum(eigen((STilde-TTilde)%*%solve(W)%*%(STilde-TTilde)%*%solve(W))$values))
    if(is.na(delta2)){stop("Test statistic can't be computed")}
    return(Re(delta2))
  }else{
    if(det(W)==0){
      return(NA)
    }else{
      delta2 <- sum(eigen((STilde-TTilde)%*%solve(W)%*%(STilde-TTilde)%*%solve(W))$values)
      return(Re(delta2))
    }
  }
}
