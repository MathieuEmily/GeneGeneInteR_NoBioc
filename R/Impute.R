#'@title Missing values handling for SnpMatrix object
#'
#'@description \code{imputeSnpMatrix} is a generic wrapper of
#'  \code{snp.imputation} and \code{impute.snps} functions from \code{snpStats}
#'  package. This function mimics a Leave-One-Out process where missing SNP are
#'  imputed for an individual based on a model trained on all other individuals.
#'  It is always preferable to use \code{snp.imputation} and \code{impute.snps}
#'  with settings specific to imputed dataset.
#'
#'@details For each individual in the dataset, following process is performed:
#'  \itemize{\item missing SNP are detected for individual i \item rules are
#'  imputed for every missing SNP using the whole dataset where individual i is
#'  removed \item SNP are imputed for individual i} This allows to impute a
#'  large number of missing SNP but some missing values may remains.
#'
#'  In that case the choice is left to the user about what action to take: leave
#'  the dataset as it is, remove all SNP or individuals with missing values.
#'  Removing all SNP is often more parsimonious than removing individuals and
#'  allows to get a dataset without any missing values with minimum
#'  information-loss.
#'
#'@param snpX \code{snpMatrix} object of which SNP are to be removed
#'@param genes.info data frame with four columns named \code{Genenames},
#'  \code{SNPnames}, \code{Position} and \code{Chromosome}. Each row describes a
#'  SNP and missing values are not allowed. If \code{genes.info} is provided an
#'  updated version is returned.
#'@param on.rem \emph{(optional)} a character string matching one of the
#'  following items: SNP, ind, none. Describes the action taken in case of
#'  remaining missing values. \emph{See details}
#'@param quiet \emph{(optional)} a boolean describing wether or not progress bar
#'  should be displayed.
#'
#'@return a list object with two named elements:\describe{ \item{\code{snpX}}{a
#'  \code{SnpMatrix} object corresponding to input matrix with imputed values.}
#'  \item{\code{genes.info}}{a data frame object corresponding to an updated
#'  version of input \code{genes.info} in case SNP had to be removed.}}
#'
#'  A warning is issued when SNP or individuals had to be removed.
#'
#'@seealso \code{\link{GGI}} \code{\link{snpMatrixScour}}
#'
#'@export
#'
#'@examples
#' ped <- system.file("extdata/example.ped", package="GGItest")
#' info <- system.file("extdata/example.info", package="GGItest")
#' posi <- system.file("extdata/example.txt", package="GGItest")
#' dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\t")
#'
#' ## In this example, genes are loosely scoured but default are much harsher
#' ## conditions
#' imputed.snps <- imputeSnpMatrix(dta$snpX, genes.info = dta$genes.info)
imputeSnpMatrix <- function(snpX, genes.info,
                            on.rem = c("SNP", "ind", "none"), quiet=FALSE){
  if (class(snpX) != "SnpMatrix") {
    stop("snpX argument should be SnpMatrix object.")
  } else if (!is.null(genes.info) && (!(is.data.frame(genes.info) | nrow(genes.info) > ncol(snpX)))) {
    stop("genes.info should be a data.frame with less rows than or as much as snpX columns.")
  } else if (!is.null(genes.info) && ncol(genes.info) != 4) {
    stop("genes.info argument should have four columns.")
  } else if (!is.null(genes.info) && !all(names(genes.info) %in% c("Genenames", "SNPnames", "Position", "Chromosome"))) {
    stop("genes.info argument should have its columns named: Genenames, SNPnames, Position, Chromosome")
  } else if (!is.null(genes.info) && is.character(genes.info$Genenames)) {
    stop("gene.info argument's Gene.name column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$SNPnames)) {
    stop("gene.info argument's SNP.name column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$Position)) {
    stop("gene.info argument's Position column should be of class character.")
  } else if (!is.null(genes.info) && is.character(genes.info$Chromosome)) {
    stop("gene.info argument's Chromosome column should be of class character.")
  } else if (!is.null(genes.info) && any(is.na(genes.info))) {
    stop("genes.info can't have missing values (NA).")
  } else if (!is.logical(quiet)) {
    stop("quiet argument should be a logical.")
  }

  on.rem <- match.arg(on.rem)

  imputed <- as(snpX, "numeric")
  if (!quiet){prog <- txtProgressBar(0, nrow(snpX), char="-", style = 3)}
  for (i in 1:nrow(snpX)) {

    select <- which(is.na(snpX[i, ]))

    if (length(select) > 0) {
      missing <- snpX[-i, select]
      present <- snpX[-i, -select]

      pos.miss <- genes.info$Position[select]
      pos.pres <- genes.info$Position[-select]

      rules <- silence(snpStats::snp.imputation)(present, missing, pos.pres, pos.miss)

      imp.targ <- snpStats::impute.snps(rules, snpX[i, ])

      imputed[i, is.na(imputed[i, ])] <- round(imp.targ)
    }

    if (!quiet){setTxtProgressBar(prog, i)}
  }
  cat("\n")

  if (any(is.na(imputed)) && on.rem == "SNP") {
    select <- which(colSums(is.na(imputed)) > 0)
    genes.info <- droplevels(genes.info[!(genes.info$SNPnames %in% colnames(imputed)[select]), ])
    imputed <- imputed[, -select]

    warning(paste(length(select), "SNP were removed due to remaining missing values."))
  } else if (any(is.na(imputed)) && on.rem == "ind") {
    select <- which(colSums(is.na(imputed)) > 0)
    imputed <- imputed[-select,]

    warning(paste(length(select), "individuals were removed due to remaining missing values."))
  } else if (any(is.na(imputed))) {
    warning(paste(sum(is.na(imputed)), "remaining missing values."))
  }

  imputed <- as(imputed, "SnpMatrix")

  return(list(snpX = imputed, genes.info = genes.info))
}

# Function used to silence snpStats::snp.imputation
silence <- function(f){
  return(function(...) {capture.output(w<-f(...));return(w);});
}
