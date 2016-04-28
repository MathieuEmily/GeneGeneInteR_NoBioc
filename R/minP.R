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
