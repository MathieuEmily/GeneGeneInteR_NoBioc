#'@title Gene-Gene Interaction Analysis Wrapper
#'
#'@description \code{GGI} performs a Gene-Gene Interaction analysis. Several
#'  methods are available among which:\itemize{\item PCA-based analysis\item
#'  CCA-based analyis\item KCCA-based analysis\item CLD-based analysis\item
#'  PLSPM-based analysis\item GBIGM \item minP test \item GATES test \item tTS
#'  test \item tProd test}
#'
#'@details This function is a wrapper for all GGI analysis methods and drive the
#'  overall analysis: splitting the dataset in genes matrices and starting
#'  elementary analysis for each pair of genes.
#'
#'  \emph{SNP from the same gene must be contiguous. See
#'  \code{\link{select.snps}}}
#'
#'  If \code{genes.lenghts} is provided, the values in it are the number of SNP
#'  from the last gene. Meaning that if \code{genes.length} is as follow: c(20,
#'  35, 15), then gene 1 will be considered to be the first 20 columns/SNP of
#'  \code{snpX}, gene 2 will be considered to be the following 35 columns/SNP,
#'  and so on...
#'
#'  If \code{genes.info} is provided instead,  then it is only important that
#'  SNP of a same gene are contiguous. Corresponding columns are automatically
#'  detected.
#'
#'  Two work scales are available of which there are specific methods.
#'
#'  \subsection{Gene's scale}{ When working at gene's scale the available
#'  methods are:\itemize{ \item Principal Components Analysis method (PCA) - PCA
#'  is performed on both genes' matrices and resulting principal components are
#'  used to fit a logit model and test genes' interaction. \item Canonical
#'  Correlation Analysis (CCA) - The maximum of canonical correlation between
#'  the two genes is computed for each group (case and control). The difference
#'  between the two transformed values (Fisher transformation) is used to test
#'  genes' interaction. \item Kernel Canonical Correlation Analysis (KCCA) -
#'  This method is the same as CCA but the canonical correlations are computed
#'  using Kernel method. \item Complete Linkage Desequilibrium (CLD) - Analysis
#'  based on the difference of the covariance matrices for case and control,
#'  between the two genes. A method based on Nagao normalized Quadratic Distance
#'  is used to compute the test statistic. \item Partial Least Square Path
#'  Modeling (PLSPM) - A network of statistical relations between latent and
#'  manifest variables is built. The difference between the path coefficients is
#'  used to compute the test statistic. \item Gene-Based Information Gain Method
#'  (GBIGM) - Entropies are used to calculate the test statistic and test genes'
#'  interaction.} }
#'
#'  \subsection{SNP's scale}{ When working at SNP's scale, a SSI is first
#'  performed beween each genes combinations and p-values are then calculated
#'  for GGI using following methods: \itemize{\item Minimum p-value test (minP)
#'  - The strongest signal is evaluated using multivariate normal distribution.
#'  \item Gene Association Test using Extended Simes procedure (GATES) - Signals
#'  are transformed and the strongest one is evaluated. \item Truncated Tail
#'  Strength test (tTS) - All signals lesser than a set threshold are
#'  transformed and evaluated. \item Truncated p-value Product test (tProd) -
#'  Same than tTS with a different signal transformation.} }
#'
#'  Missing values are rejected and trying to parse an incomplete
#'  \code{SnpMatrix} object as an argument will result in an error.
#'
#'@param Y numeric or factor vector with two values (most often 0, 1). This is
#'  the response variable and should be of length equal to the number of row of
#'  \code{G1} and \code{G2} (number of individuals).
#'@param snpX \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix}
#'  object. Must have a number of row equal to Y argument's length.
#'@param genes.length numeric vector. It is the length (in columns/SNP) of each
#'  gene. Each gene declared is considered contiguous with the one before and
#'  after it. \code{genes.length} can be named if you want the returned matrix
#'  to have dimensions named after those. If no names are given then generic
#'  names are generated following the pattern Gene.n (n being the gene's index)
#'  .
#'@param genes.info \emph{(optional)} a data frame. It must have four columns
#'  named \code{Genenames}, \code{SNPnames}, \code{Position} and
#'  \code{Chromosome}. Each row describes a SNP and missing values are not
#'  allowed.
#'@param method a string matching onz of the following: PCA.Std, PCA.GenFreq,
#'  CCA, KCCA, CLD, PLSPM, GBIGM, minP, GATES, tTS or tProd. Only one method can
#'  be parsed.
#'@param ... Other optional arguments to be passed to the functions associated
#'  with the method chosen. See more in elementary methods help.
#'
#'@return \code{GGI} returns a symmetric matrix of size \eqn{G*G} where \eqn{G}
#'  is the number of genes studied. The general term of the matrix is the
#'  p-value of the interaction between the two genes.
#'
#'@seealso \code{\link{PCA.Std}}, \code{\link{PCA.GenFreq}},
#'  \code{\link{CCA.test}}, \code{\link{KCCA.test}}, \code{\link{CLD.test}},
#'  \code{\link{PLSPM.test}}, \code{\link{GBIGM.test}}, \code{\link{GGI.plot}},
#'  \code{\link{minP.test}}, \code{\link{gates.test}}, \code{\link{tTS.test}},
#'  \code{\link{tProd.test}}
#'
#'@export
#'
#'@examples
#' ## Dataset is included in the package
#' ped <- system.file("extdata/example.ped", package="GGItest")
#' info <- system.file("extdata/example.info", package="GGItest")
#' posi <- system.file("extdata/example.txt", package="GGItest")
#'
#' dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\t")
#' dta <- imputeSnpMatrix(dta$snpX, genes.info = dta$genes.info)
#' resp <- system.file("extdata/response.txt", package="GGItest")
#'
#' # If a data frame is provided to GGI or one of the *.test function, only the
#' # column is checked and used.
#' Y  <- read.csv(resp, header=FALSE)
#'
#' ## By default the PCA-based method is used.
#' GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info)
GGI <- function(Y, snpX, genes.length = NULL, genes.info = NULL,
            method = c("PCA.Std", "PCA.GenFreq", "CCA", "KCCA","CLD","PLSPM","GBIGM",
                       "minP", "GATES", "tTS", "tProd"), ...){

  if (!is.null(dim(Y))) {
    Y <- Y[, 1]
  }

  # Arguments checks
  if (nlevels(as.factor(Y)) != 2) {
    stop("response variable should be binary. (2 modes).")
  } else if (class(snpX) != "SnpMatrix") {
    stop("snpX argument should be SnpMatrix object.")
  } else if (length(Y) != nrow(snpX)) {
    stop("Response variable should be conformant with genes matrix rows number.")
  } else if (is.null(genes.length) && is.null(genes.info)) {
    stop("Genes must be defined with either genes.length or genes.info arguments.")
  } else if (!is.null(genes.length) && !is.numeric(genes.length)) {
    stop("genes.length should be a numeric vector.")
  } else if (!is.null(genes.length) && any(genes.length <= 0)) {
    stop("Genes length can't be lesser or equal to 0.")
  } else if (!is.null(genes.length) && length(genes.length) < 2) {
    stop("Select at least two genes.")
  } else if (!is.null(genes.length) &&  sum(genes.length) > ncol(snpX)) {
    stop("genes.length argument non conformant with snps.matrix rows count.")
  } else if (!is.null(genes.length) && any(is.na(genes.length))) {
    stop("When provided, genes.length argument can't have missing values (NA).")
  } else if (!is.null(genes.info) && (!(is.data.frame(genes.info) | nrow(genes.info) > ncol(snpX)))) {
    stop("When provided, genes.info should be a data.frame with less rows than or as much as snpX columns.")
  } else if (!is.null(genes.info) && ncol(genes.info) != 4) {
    stop("genes.info argument should have four columns. See help file.")
  } else if (!is.null(genes.info) && !all(names(genes.info) %in% c("Genenames", "SNPnames", "Position", "Chromosome"))) {
    stop("genes.info argument should have its columns named: Genenames, SNPnames, Position, Chromosome")
  } else if (!is.null(genes.info) && is.character(genes.info$Genenames)) {
    stop("gene.info argument's Gene.name column should be of class character.")
  } else if (!is.null(genes.info) && nlevels(genes.info$Genenames) < 2) {
    stop("Select at least two genes.")
  } else if (!is.null(genes.info) && is.character(genes.info$SNPnames)) {
    stop("gene.info argument's SNP.name column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$Position)) {
    stop("gene.info argument's Position column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$Chromosome)) {
    stop("gene.info argument's Chromosome column should be of class character.")
  } else if (!is.null(genes.info) && any(is.na(genes.info))) {
    stop("When provided, genes.info can't have missing values (NA).")
  } else if (any(is.na(snpX))) {
    stop("The snpMatrix must be complete. No NAs are allowed.")
  } else if (any(is.na(Y))) {
    stop("The response variable must be complete. No NAs are allowed.")
  } else if (!is.null(genes.length) && !is.null(genes.info)) {
    warning("Both genes.length and genes.info arguments were provided, only genes.info will be used.")
  }

  method <- match.arg(method)

  Y <- as.numeric(Y)
  if(min(Y)!=0){Y<-Y-min(Y)}
  Y <- as.factor(Y)

  #SnpMatrix and dataframe are reordered to make sure that all SNP of a gene are contiguous.
  if (!is.null(genes.info)){
    genes.info <- genes.info[order(genes.info$Genenames, genes.info$SNPnames), ]
    snpX <- snpX[, as.character(genes.info$SNPnames)]
  }

  #Genes names computing
  if (!is.null(genes.info)) {
    genes.names <- levels(genes.info$Genenames)
  } else if (!is.null(genes.length) && is.null(names(genes.length))){
    genes.names <- paste("Gene", seq(1:length(genes.length)), sep=".")
    names(genes.length) <- genes.names
  } else {
    genes.names <- names(genes.length)
  }

  #Interactions listing
  interactions <- combn(genes.names, m=2)

  #Indexes of the genes
  if (!is.null(genes.info)) {
    gene.start <- NULL
    gene.end <- NULL
    for (i in 1:nlevels(genes.info$Genenames)) {
      gene <- genes.info[which(genes.info$Genenames %in% levels(genes.info$Genenames)[i]), ]
      gene.start <- c(gene.start, min(which(colnames(snpX) %in% gene$SNPnames), na.rm = TRUE))
      gene.end   <- c(gene.end, max(which(colnames(snpX) %in% gene$SNPnames), na.rm = TRUE))
    }
  } else {
    gene.start <- c(0, cumsum(genes.length)[1:(length(genes.length) - 1)]) + 1
    gene.end   <- cumsum(genes.length)
  }

  names(gene.start) <- genes.names
  names(gene.end)   <- genes.names

  #Setup of the return object
  genes.interactions <- diag(0, length(genes.names))
  colnames(genes.interactions) <- genes.names
  rownames(genes.interactions) <- genes.names

  #Application of the method on the interactions
  for (i in 1:ncol(interactions)) {
    print(paste("Interaction between", interactions[1, i], "&", interactions[2, i]))

    G1 <- snpX[, gene.start[interactions[1, i]]:gene.end[interactions[1, i]]]
    G2 <- snpX[, gene.start[interactions[2, i]]:gene.end[interactions[2, i]]]

    if(!method %in% c("minP","GATES","tTS","tProd") || ncol(G1)*ncol(G2)<1000){
    	print(i)
    genes.interactions[interactions[1, i], interactions[2, i]] <- switch(method,
                                                                         CCA = CCA.test(Y, G1, G2, ...)$p.value,
                                                                         KCCA = KCCA.test(Y, G1, G2, ...)$p.value,
                                                                         CLD = CLD.test(Y, G1, G2, ...)$p.value,
                                                                         PLSPM = PLSPM.test(Y, G1, G2, ...)$p.value,
                                                                         GBIGM = GBIGM.test(Y, G1, G2, ...)$p.value,
                                                                         PCA.Std = PCA.Std(Y, G1, G2, ...)$p.value,
                                                                         PCA.GenFreq = PCA.GenFreq(Y, G1, G2, ...)$p.value,
                                                                         minP = minP.test(Y, G1, G2, ...)$p.value,
                                                                         GATES = gates.test(Y, G1, G2, ...)$p.value,
                                                                         tTS   = tTS.test(Y, G1, G2, ...)$p.value,
                                                                         tProd = tProd.test(Y, G1, G2, ...)$p.value)
    } else {
      warning("Too much interactions to test for SSI method (>1000). NA returned")
      genes.interactions[interactions[1, i], interactions[2, i]] <- NA
    }
  }

  genes.interactions[lower.tri(genes.interactions)] <- t(genes.interactions)[lower.tri(genes.interactions)]
  return(genes.interactions)
}

# Function that order genes.info and SnpMatrix so that all SNP of a genes are contiguous
order.snpMatrix <- function(snpX, genes.info) {
  genes.info <- genes.info[order(genes.info$Genenames, genes.info$SNPnames), ]
  snpX <- snpX[, as.character(genes.info$SNPnames)]

  return(list(snpX = snpX, genes.info = genes.info))
}
