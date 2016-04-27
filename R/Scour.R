#'@title SNP filtering based on MAF and Hardy-Weinberg Equilibrium
#'
#'@description \code{snpMatrixScour} filters SNP of a \code{snpMatrix} object
#'  based on Minor Allele Frequency criterion and Hardy-Weinberg Equilibrium.
#'  Filtering criteria can be adjusted by the user.
#'
#'@details This function removes SNP that are not conformant with Minor Allele
#'  Frequency criterion (MAF) or Hardy-Weinberg equilibrium (HWE) criterion.
#'  \itemize{ \item Every SNP whose MAF is inferior to user-set threshold is
#'  removed. Default value is 1\%. \item Every SNP that does not verify HWE is
#'  removed. Deviation to HWE is calculated with a \eqn{\chi ^2} test and SNP is
#'  removed if resulting p-value is less than user-set threshold. Default value
#'  is 1\%.}
#'
#'  If gene lengths are provided, as a numeric vector or as a data frame, they
#'  are updated to keep track of SNP removing. Those objects can then be
#'  directly used with \code{GGI} function.
#'
#'  Missing values are rejected and trying to parse an incomplete
#'  \code{SnpMatrix} object as an argument will result in an error.
#'
#'@param snpX \code{snpMatrix} object of which SNP are to be removed
#'@param genes.length \emph{(optional)} numeric vector. It is the length (in
#'  columns/SNP) of each gene. Each gene declared is considered contiguous with
#'  the one before and after it. \code{genes.lengths} can be named (names will
#'  be kept). If \code{genes.length} is provided an updated version is returned.
#'@param genes.info \emph{(optional)} a data frame. It must have four columns
#'  named \code{Genenames}, \code{SNPnames}, \code{Position} and
#'  \code{Chromosome}. Each row describes a SNP and missing values are not
#'  allowed. If \code{genes.info} is provided an updated version is returned.
#'@param min.maf a numeric that is the minimum MAF (Minor Allele Frequency) for
#'  a SNP. SNP that does not meet that criterion are removed. Default is 1\%.
#'@param min.eq a numeric that is the maximum acceptable p-value for the
#'  \eqn{\chi ^2} verifying HWE deviation. SNP that does not meet that criterion
#'  are removed. Default is 1\%.
#'@param NA.rate a numeric that is the maximum acceptable frequency of NA values
#'  . Default is 10\%. High frequencies of missing values (NA) can make
#'  imputation harder (residual missing values).
#'
#'@return A \code{SnpMatrix} object is always returned. If gene lengths were
#'  provided then an updated object of the same class is also returned, in that
#'  case both the SnpMatrix and the gene lengths object are returned in a named
#'  list:\describe{\item{\code{snpX}}{the \code{SnpMatrix} of which
#'  non-conformant SNP were removed.} \item{\code{genes}}{the object that
#'  contains gene lengths information. Can be a numeric vector (possibly named)
#'  or a data frame.}}
#'
#'  A warning message is issued when a list is returned.
#'
#'@seealso \code{\link{GGI}}
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
#' new.snps <- snpMatrixScour(dta$snpX, genes.info = dta$genes.info,
#'                            min.maf = 0.2, min.eq=0.05, NA.rate = 0.2)

snpMatrixScour <- function(snpX, genes.length = NULL, genes.info = NULL,
                           min.maf = 0.01, min.eq = 0.01, NA.rate=0.1) {
  if (class(snpX) != "SnpMatrix") {
    stop("snpX argument should be SnpMatrix object.")
  } else if (!is.null(genes.length) && (!(is.numeric(genes.length) || sum(genes.length) > ncol(snpX)))) {
    stop("genes.length argument should be a numeric vector which sum should be lesser than ncol(snpX).")
  } else if (!is.null(genes.length) && any(is.na(genes.length))) {
    stop("genes.length argument can't have missing values (NA).")
  } else if (!is.null(genes.info) && (!(is.data.frame(genes.info) || nrow(genes.info) > ncol(snpX)))) {
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
  } else if (!is.numeric(min.maf) || min.maf < 0 || min.maf > 1) {
    stop("min.maf argument should be a numeric between 0 & 1.")
  } else if (!is.numeric(min.eq) || min.eq < 0 || min.eq > 1) {
    stop("min.eq argment should be a numeric between 0 & 1.")
  } else if (!is.numeric(NA.rate) || NA.rate > 1 || NA.rate < 0) {
    stop("NA.rate argment should be a numeric between 0 & 1.")
  } else if (!is.null(genes.length) && !is.null(genes.info)) {
    stop("Both genes.length and genes.info arguments were provided, only genes.info will be used.")
  }

  # If no genes were defined, whole dataset is scoured as is
  if (is.null(genes.length) && is.null(genes.info)) {
    return(GeneScour(snpX, min.maf, min.eq, NA.rate))
  # If genes are defined, each genes is scoured individually to
  # keep track of its length variation.
  } else {
    #Indexes of the genes
    if (!is.null(genes.info)) {

      #SnpMatrix and dataframe are reordered to make sure that all SNP of a gene are contiguous.
      genes.info <- genes.info[order(genes.info$Genenames, genes.info$SNPnames), ]
      snpX <- snpX[, as.character(genes.info$SNPnames)]

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

    n.genes <- length(gene.start)
    new.snpX <- NULL
    for (i in 1:n.genes) {
      cur.gene <- try(GeneScour(snpX[, gene.start[i]:gene.end[i]], min.maf, min.eq, NA.rate), silent=TRUE)

      if (class(cur.gene) == "try-error") {
        warning("A whole gene had to be removed as no SNP met the requirements.")
      }

      # Gene positions are updated
      if (is.null(genes.info)) {
        if (class(cur.gene) == "try-error"){
          new.genes[i] <- NA
        } else {
          new.genes[i] <- ncol(cur.gene)
        }
      } else {
        if (class(cur.gene) == "try-error"){
          cur.gene.name <- levels(genes.info$Genenames)[i]
          genes.info <- genes.info[- which(genes.info$Genenames == cur.gene.name),]
        } else {
          SNP.cond  <- !(genes.info$SNPnames %in% colnames(cur.gene))
          gene.cond <- genes.info$Genenames == levels(genes.info$Genenames)[i]
          select <- which(SNP.cond & gene.cond)
          if (length(select) > 0) {genes.info <- genes.info[- select,]}
        }
      }

      # Binding matrices
      if (class(cur.gene) != "try-error" && is.null(new.snpX)) {
        new.snpX <- cur.gene
      } else if (class(cur.gene) != "try-error") {
        new.snpX <- cbind(new.snpX, cur.gene)
      }
    }

    if (!is.null(genes.info)) {
      new.genes <- droplevels(genes.info)
    } else if (length(which(is.na(new.genes))) > 0) {
      new.genes <- new.genes[-which(is.na(new.genes))]
    }
  }

  new.snpX <- as(new.snpX, "SnpMatrix")

  print("A list object has been returned with elements: snpX & genes")
  return(list(snpX = new.snpX, genes = new.genes))
}

# Function that checks SNP validity
# A SNP is considered valid if it meets following criteria:
#     - MAF is superior to min.maf arguments
#     - Hardy-Weinberg equilibrium is verified using min.eq threshold
GeneScour <- function(gene, min.maf = 0.01, min.eq = 0.01, NA.rate=0.1){
  # NA filter
  SNP.NA <- colSums(is.na(gene))/nrow(gene)

  # MAF filtering
  SNP.MAF <- snpStats::col.summary(gene)$MAF

  # Hardy-Weinberg Equilibrium check
  AA <- colSums(as(gene, "numeric") == 2, na.rm=TRUE)
  Aa <- colSums(as(gene, "numeric") == 1, na.rm=TRUE)
  aa <- colSums(as(gene, "numeric") == 0, na.rm=TRUE)

  p <- (2*AA + Aa)/(2*(AA + Aa + aa))
  q <- 1 - p

  E.AA <- p^2 * nrow(gene)
  E.Aa <- 2 * p * q * nrow(gene)
  E.aa <- q^2 * nrow(gene)

  # Chi-Square test
  X.AA <- ((AA - E.AA)^2)/E.AA
  X.Aa <- ((Aa - E.Aa)^2)/E.Aa
  X.aa <- ((aa - E.aa)^2)/E.aa
  X <- X.AA + X.Aa + X.aa
  p.X <- pchisq(X, 1, lower.tail=FALSE)

  # Filtering
  if (any(SNP.MAF >= min.maf & p.X >= min.eq & SNP.NA <= NA.rate)) {
    return(gene[, (SNP.MAF >= min.maf & p.X >= min.eq & SNP.NA <= NA.rate)])
  } else {
    stop("No SNP with sufficient MAF and verifying Hardy-Weinberg equilibrium")
  }
}
