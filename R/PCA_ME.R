#'@title PCA-based GGI analysis.
#'
#'@description \code{PCA.Std} performs a Gene-Gene Interaction (GGI) analysis
#'  using a PCA in which dataset is standardized using variables' standard
#'  deviation. With \code{PCA.GenFreq}, dataset is standardized using the
#'  standard deviation under Hardy-Weinberg equilibrium.
#'
#'@details A Principal Components Analysis is performed on both \code{G1} and
#'  \code{G2} matrices. With \code{PCA.Std} dataset is standardized using
#'  variables' standard deviation. Whereas with \code{PCA.GenFreq}, dataset is
#'  standardized using standard deviation under Hardy-Weinberg equilibrium. Each
#'  of these two matrices represents a gene to test and is built as follow:
#'  \itemize{ \item A column is a SNP of said gene. \item A row is an individual
#'  for which the haplotype and the response variable are known.}
#'
#'  Principal components are retrieved to describe each dataset with user-set
#'  inertia percentage.
#'
#'  Principal components are then used in a logit model fitting process in which
#'  interaction is tested. If model fitting fails then a PC is removed from the
#'  pool and a new attempt to fitting is made. The process loops until either
#'  the model is fitted or there are only two PC for each genes remaining.
#'  Principal components are removed according to lowest inertia percentage.
#'
#'@param Y numeric or factor vector with two values (most often 0, 1). This is
#'  the response variable and should be of length equal to the number of row of
#'  \code{G1} and \code{G2} (number of individuals).
#'@param G1 \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
#'  Must have a number of row equal to Y argument's length.
#'@param G2 \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object.
#'  Must have a number of row equal to Y argument's length.
#'@param threshold (optional) numeric comprised in \eqn{[0, 1]} or
#'  \eqn{[0, 100]} intervals. This is how much inertia should be kept for each
#'  PCA and has a direct impact on how many PCA dimensions are kept.
#'
#'@return If the function succeed, the interaction p-value is returned.
#'
#'  A warning can be issued if the inertia threshold is too high and principal
#'  components had to be removed.
#'
#'  An error is returned if too much PCA dimensions had to be removed.
#'
#'@references Jia Li, Rui Tang, Joanna Biernacka and Maricka de Andrade.
#'  (2009) Identification of gene-gene interaction using principal components. BMC
#'  Proceedings, 3 (Suppl. 7): S78
#'
#'@seealso \code{\link{GGI}}
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
#' dta <- imputeSnpMatrix(dta$snpX, genes.info = dta$genes.info)
#'
#' G1 <- select.snps(dta$snpX, dta$genes.info, "TXNDC5")$snpX
#' G2 <- select.snps(dta$snpX, dta$genes.info, "DNAH9")$snpX
#'
#' resp <- system.file("extdata/response.txt", package="GGItest")
#' Y  <- read.csv(resp, header=FALSE)
#'
#' ## By default, inertia threshold is 80%
#' PCA.Std(Y, G1, G2)
#'
#' ## Setting a higher inertia threshold often ends up with removed PCA
#' ## dimensions
#' PCA.GenFreq(Y, G1, G2, threshold = 1)
PCA.Std <- function(Y, G1, G2, threshold=0.8) {

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  # Arguments checks
  if (class(threshold) != "numeric") {
    stop("thresold argument should be a numeric.")
  } else if (threshold < 0 | threshold > 100) {
    stop("threshold argument shoud be comprised in [0, 1] or [0, 100] interval.")
  } else if (nlevels(as.factor(Y)) != 2) {
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

  # Threshold formatting
  if (threshold < 1) {
    inertia.thresh <- threshold * 100
  } else {
    inertia.thresh <- threshold
  }

  # SnpMatrix coerced into matrix to be compatible with FactoMineR::PCA
  G1.num <- as(G1, "numeric")
  G2.num <- as(G2, "numeric")

  # PCA performed
  G1.PCA <- FactoMineR::PCA(G1.num, ncp=NULL, graph=FALSE)
  G2.PCA <- FactoMineR::PCA(G2.num, ncp=NULL, graph=FALSE)

  # Genes are represented by PCA coords on enough dimensions to retrieve
  # as much inertia as set by user.
  G1.ncp <- which(G1.PCA$eig[, 3] > inertia.thresh)[1]
  G2.ncp <- which(G2.PCA$eig[, 3] > inertia.thresh)[1]

  if (G1.ncp == 1) {
    G1.VarSynth <- data.frame(Dim.1 =G1.PCA$ind$coord[, 1:G1.ncp])
  } else {
    G1.VarSynth <- G1.PCA$ind$coord[, 1:G1.ncp]
  }

  if (G2.ncp == 1) {
    G2.VarSynth <- data.frame(Dim.1 =G2.PCA$ind$coord[, 1:G2.ncp])
  } else {
    G2.VarSynth <- G2.PCA$ind$coord[, 1:G2.ncp]
  }

  G1.PCA <- list(VarSynth = G1.VarSynth, Inertia = G1.PCA$eig[1:G1.ncp, 2])
  G2.PCA <- list(VarSynth = G2.VarSynth, Inertia = G2.PCA$eig[1:G2.ncp, 2])

  # Interaction effects are tested
  return(compare.PCA(Y, G1.PCA, G2.PCA))
}

#'@describeIn PCA.Std Standardization based on Hardy-Weinberg equilirum
#'@export
PCA.GenFreq <- function(Y, G1, G2, threshold=0.8) {

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  # Arguments checks
  if (class(threshold) != "numeric") {
    stop("thresold argument should be a numeric.")
  } else if (threshold < 0 | threshold > 100) {
    stop("threshold argument shoud be comprised in [0, 1] or [0, 100] interval.")
  } else if (nlevels(as.factor(Y)) != 2) {
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

  # Threshold formatting
  if (threshold < 1) {
    inertia.thresh <- threshold * 100
  } else {
    inertia.thresh <- threshold
  }

  G1.PCA <- get.PCA.res(G1)
  G2.PCA <- get.PCA.res(G2)

  # Genes are represented by PCA coords on enough dimensions to retrieve
  # as much inertia as set by user.
  G1.ncp <- which(cumsum(G1.PCA$Inertia) > inertia.thresh)[1]
  G2.ncp <- which(cumsum(G2.PCA$Inertia) > inertia.thresh)[1]

  G1.PCA <- list(VarSynth=G1.PCA$VarSynth[, 1:G1.ncp], Inertia=G1.PCA$Inertia[1:G1.ncp])
  G2.PCA <- list(VarSynth=G2.PCA$VarSynth[, 1:G2.ncp], Inertia=G2.PCA$Inertia[1:G2.ncp])

  if (G1.ncp == 1) {
    G1.PCA$VarSynth <- data.frame(Dim.1 = G1.PCA$VarSynth)
  }

  if (G2.ncp == 1) {
    G2.PCA$VarSynth <- data.frame(Dim.1 = G2.PCA$VarSynth)
  }

  # Interaction effects are tested
  return(compare.PCA(Y, G1.PCA, G2.PCA))
}

## Function that retrieves Principal Components and Eigen Values
## of a SnpMatrix object.
get.PCA.res <- function(gene.matrix){
  if (class(gene.matrix) != "SnpMatrix") {
    stop("gene.matrix argument should be SnpMatrix object.")
  } else if (sum(is.na(gene.matrix))!=0) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  }

  Id <- diag(ncol(gene.matrix))
  # X'X matrix is calculated
  G.XpX <- snpStats::snp.pre.multiply(gene.matrix, t(snpStats::snp.post.multiply(gene.matrix, Id)))
  # Singular Values Decomposition is performed
  G.eigen <- eigen(G.XpX, symmetric=TRUE)
  # Individuals coordinates are calculated
  G.PCA <- snpStats::snp.post.multiply(gene.matrix, G.eigen$vectors)
  # Columns names are formatted (for comparison procedure)
  colnames(G.PCA) <- paste("Dim", seq(1, ncol(G.PCA)), sep=".")

  # Computes inertia percentages for all ordered eigen values
  eigen.val <- 100*G.eigen$values/sum(G.eigen$values)

  return(list(VarSynth=G.PCA, Inertia=eigen.val))
}

## Function that fits single effects model and interaction
## model and compare them.
compare.PCA <- function(Resp., G1.PCA, G2.PCA) {
  # Arguments checks
  if (nlevels(as.factor(Resp.)) != 2) {
    stop("response variable should be binary. (2 modes).")
  } else if (class(G1.PCA) != "list" | class(G2.PCA) != "list") {
    stop("G1.PCA and G2.PCA arguments should be list objects.")
  } else if (length(G1.PCA) != 2 | length(G2.PCA) != 2) {
    stop("G1.PCA and G2.PCA arguments should be of length 2: VarSynth and Inertia elements.")
  } else if (nrow(G1.PCA$VarSynth) != nrow(G2.PCA$VarSynth)) {
    stop("VarSynth elements of G1.PCA and G2.PCA should have same rows count.")
  } else if (length(Resp.) != nrow(G1.PCA$VarSynth) | length(Resp.) != nrow(G2.PCA$VarSynth)) {
    stop("Response variable should be conformant with genes matrices rows number.")
  } else if (length(G1.PCA$Inertia) > ncol(G1.PCA$VarSynth) | length(G2.PCA$Inertia) > ncol(G2.PCA$VarSynth)){
    stop("Inertia elements of G1.PCA and G2.PCA can't be of greater length than ncol(VarSynth)")
  } else if (sum(is.na(Resp.)) != 0) {
    stop("The response variable vector must be complete. No NAs are allowed.")
  }

  G1.inertia <- G1.PCA$Inertia
  G2.inertia <- G2.PCA$Inertia
  G1.PCA <- G1.PCA$VarSynth
  G2.PCA <- G2.PCA$VarSynth

  # Effects' names
  colnames(G1.PCA) <- paste("G1", colnames(G1.PCA), sep=".")
  colnames(G2.PCA) <- paste("G2", colnames(G2.PCA), sep=".")

  # Iterating until model is fitted
  # Less informant component is removed each time fitting is impossible.
  keep.on <- TRUE
  trimmed <- FALSE
  while (keep.on) {
    data <- data.frame(Resp., G1.PCA, G2.PCA)

    # Effects set up
    single.effects <- paste(colnames(G1.PCA), colnames(G2.PCA), sep="+", collapse="+")

    # Creating all possible pairs
    inter.effects <- expand.grid(colnames(G1.PCA), colnames(G2.PCA))
    # Creating a string with all interactions
    inter.effects <- paste(inter.effects[, 1], inter.effects[, 2], sep=":", collapse="+")

    inter.formula <- as.formula(paste("Resp.~", single.effects, "+", inter.effects, sep=""))
    null.formula <- as.formula(paste("Resp.", single.effects, sep="~"))

    # Trying to fit the largest model
    inter.mod <- tryCatch(glm(inter.formula, data=data, family = "binomial"),
             warning = function(w){"warning"},
             error   = function(e){"error"})

    # If succeeded end the process
    if (any(class(inter.mod) == "glm")){
      keep.on <- FALSE
      null.mod  <- glm(null.formula, data=data, family = "binomial")
      comp.res <- anova(null.mod, inter.mod, test="Chisq")
    # If model coulnd't be fitted then a component is removed
    } else {
      trimmed <- TRUE
      # Looking for the less informant component
      if (G1.inertia[length(G1.inertia)] < G2.inertia[length(G2.inertia)]) {
		# If less informant component is one of the two first PC of the gene
		# then a PC is removed from the other gene instead.
        if (length(G1.inertia) > 2) {
          G1.inertia <- G1.inertia[-length(G1.inertia)]
          G1.PCA <- G1.PCA[, -ncol(G1.PCA)]
        } else if (length(G2.inertia) > 2){
          G2.inertia <- G2.inertia[-length(G2.inertia)]
          G2.PCA <- G2.PCA[, -ncol(G2.PCA)]
		# If both genes are only described by their first two PC then an error is issued.
        } else {
          stop("Genes too correlated to fit glm model.")
        }
      } else {
        if (length(G2.inertia) > 2) {
          G2.inertia <- G2.inertia[-length(G2.inertia)]
          G2.PCA <- G2.PCA[, -ncol(G2.PCA)]
        } else if (length(G1.inertia) > 2){
          G1.inertia <- G1.inertia[-length(G1.inertia)]
          G1.PCA <- G1.PCA[, -ncol(G1.PCA)]
        } else {
          stop("Genes too correlated to fit glm model.")
        }
      }
    }
  }

  # When PC had to be removed, a warning is issued to inform the user of
  # inertia loss.
  if (trimmed) {
    warn.str <- paste("Less principal components had to be kept to fit glm model: ",
                      round(max(cumsum(G1.inertia)),2), "% of inertia was kept for G1 & ",
                      round(max(cumsum(G2.inertia)),2), "% of inertia was kept for G2.", sep="")
    warning(warn.str)
  }

pval <- comp.res$'Pr(>Chi)'[2]
stat <- comp.res$'Deviance'[2]
df <- comp.res$'Df'[2]
	names(stat)="Deviance"
	res <- list(statistic=stat,p.value=pval,method="Principal Component Analysis",df=df)
	class(res) <- "GGItest"
  return(res)
#  return(comp.res$'Pr(>Chi)'[2])
}
