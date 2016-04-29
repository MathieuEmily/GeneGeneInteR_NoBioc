minP.test <- function(Y, G1, G2){

  # Checking if gene splitting is needed
  if (ncol(G1) * ncol(G2) > 1000){
    # Finding biggest gene and smallest gene
    if (ncol(G1) >= ncol(G2)){
      big   <- G1
      small <- G2
    } else {
      big   <- G2
      small <- G1
    }

    # Biggest gene is splitted down first, max cluster size is checked.
    # Max cluster size is the biggest cluster size that passes the filter
    # so as to keep the test number to a minimum.
    max.size <- ncol(big) - which((ncol(big):1 * ncol(small)) <= 1000)[1] + 1

    # If no no value was small enough then the smallest needs to be
    # splitted too.
    if (is.na(max.size) | max.size < 10){
      # Second gene needs to be splitted too.
      # The best combination is looked for.
      # Decreasing values of the smalled gene are tested for each of which decreasing
      # values of the smallest genes are parsed. Thus biggest gene is decreasing faster.
      max.size <- which(ncol(big):1 %o% ncol(small):1 > 1000, arrr.ind=TRUE)[1, ]

      if (is.na(max.size)){
        stop("No fitting cluster sizes could be found.")
      }

      # Converting the matrix coordinates back to cluster size
      max.size <- c(big = ncol(biggest) - max.size[1] + 1, small = ncol(smallest) - max.size[2] + 1)
    } else {
      # Splitting the biggest gene is enough
      max.size <- c(big = max.size, small = NA)
    }

    breaks.big <- find.breaks(big, max.size['big'])
    breaks.small <- find.breaks(small, max.size['small'])
  } else {
    # Arbitrarilly affecting genes to variables
    big <- G1
    small <- G2

    # Breaks are set up so that genes are parsed as one block
    breaks.big <- c(0, ncol(big))
    breaks.small <- c(0, ncol(small))
  }

  # Genes are splitted if needed
  # All pairs of sub-matrices are to be tested
  sub.pairs <- expand.grid(1:(length(breaks.big) - 1), 1:(length(breaks.small) - 1))

  # All pairs are iterated over
  pairs.p.val <- rep(NA,times=nrow(sub.pairs))
  pairs.stat <- rep(NA,times=nrow(sub.pairs))
  for (i in 1:nrow(sub.pairs)){
    # Cluster boundaries for both genes
    c.bboundaries <- (breaks.big[sub.pairs[i, 1]] + 1):breaks.big[sub.pairs[i, 1] + 1]
    c.sboundaries <- (breaks.small[sub.pairs[i, 2]] + 1):breaks.small[sub.pairs[i, 2] + 1]
    c.test <-minP.test.2pairs(Y, big[, c.bboundaries], small[, c.sboundaries])
    pairs.p.val[i] <- c.test$p.value
    pairs.stat[i] <- c.test$statistic
    
  }

	tmp <- p.adjust(pairs.p.val, "BH")
	pval <- min(tmp)[1]
	stat <- pairs.stat[which.min(tmp)]
	names(stat)="minP"
	res <- list(statistic=stat,p.value=pval,method="minP")
	class(res) <- "GGItest"
  return(res)

#  pairs.res <- min(p.adjust(pairs.res, "BH"))[1]

#  return(pairs.res)
}

find.breaks <- function(gene, clust.size){
  if (is.na(clust.size)){
    return(c(0, ncol(gene)))
  }

  # Constructing distance matrix
  distance <- snpStats::ld(gene, gene, stats='R.squared')
  distance <- as.dist(1 - distance)

  # Constructing clustering tree
  clust.tree <- rioja::chclust(distance)

  # Cutting tree
  k.clust <- cutree(clust.tree, k=ncol(gene):1)

  # Finding biggest clust size corresponding to cirteria
  for (i in ncol(gene):1){
    if (max(summary(as.factor(k.clust[, i]))) < clust.size) {
      break
    }
  }

  # At this point, i is affected with the correct value
  # summary(as.factor(k.clust[, i])) return the count for each cluster
  # As clustering was done with contiguity constraints, cumsum gives the upper bound of
  # each class. Each upper bound is the lower bound of next class, except for the first one,
  # thus 0 is added.
  breaks <- c(0, cumsum(summary(as.factor(k.clust[, i]))))

  return(breaks)
}

minP.test.2pairs <- function(Y, G1, G2){

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
