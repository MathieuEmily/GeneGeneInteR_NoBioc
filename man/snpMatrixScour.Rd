\name{snpMatrixScour}
\alias{snpMatrixScour}
\title{SNP filtering based on MAF and Hardy-Weinberg Equilibrium}
\usage{
snpMatrixScour(snpX, genes.length = NULL, genes.info = NULL,
  min.maf = 0.01, min.eq = 0.01, NA.rate = 0.1)
}
\arguments{
\item{snpX}{\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object from which SNPs are to be removed}

\item{genes.length}{\emph{(optional)} numeric vector. It is the length (in columns/SNPs) of each gene. Each gene declared is considered contiguous with the one before and after it. \code{genes.lengths} can be named (names will be kept).}

\item{genes.info}{\emph{(optional)} a \code{\link[base]{data.frame}} with four columns named \code{Genenames}, \code{SNPnames}, \code{Position} and \code{Chromosome}. Each row describes a SNP and missing values are not allowed.}

\item{min.maf}{a single value between O and O.5 that gives the threshold for the MAF (Minor Allele Frequency) of a SNP. SNP with MAF < \code{min.maf} are removed. Default is 0.01.}

\item{min.eq}{a single value between 0 and 1 that gives the maximum acceptable p-value for the \eqn{\chi ^2} verifying HWE deviation. SNP that does not meet that criterion are removed. Default is 0.01.}

\item{NA.rate}{a single value between 0 and 1 that gives the maximum acceptable proportion of missing values. Default is 0.1. High frequencies of missing values (NA) can make imputation harder (residual missing values).}
}
\value{
A list with two objects:
\describe{
\item{\code{snpX}}{the \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object where non-conform SNPs are removed.}
\item{\code{genes}}{the object that contains the updated gene lengths information. Can be a numeric vector (possibly named) or a data frame. If \code{genes.length} and \code{genes.info} are not provided by the user as input of the \code{snpMatrixScour} function, the \code{genes} object if \code{NULL}.}
}
\description{
\code{snpMatrixScour} aims at filtering out SNPs of a {\href{http://bioconductor.org/packages/snpStats/}{snpMatrix} object based on Minor Allele Frequency criterion and deviation to Hardy-Weinberg Equilibrium.
}
\details{
This function removes SNPs that does not meet all following criteria:
\itemize{
	\item \eqn{MAF > \code{min.maf}},
	\item \eqn{HWE > \code{min.eq}}, where HWE is the p-value of the \eqn{\chi ^2} test of deviation to Hardy-Weinberg Equilibrium,
	\item \eqn{NA > \code{NA.rate}}, where NA is the proportion of missing values.
}

If \code{genes.length} and \code{genes.info} are provided by the user, an updated version is returned by \code{snpMatrixScour}. The returned object can be directly used as inputs of the \code{\link{GGI}} function.
}
\examples{
ped <- system.file("extdata/example.ped", package="GeneGeneInteR")
info <- system.file("extdata/example.info", package="GeneGeneInteR")
posi <- system.file("extdata/example.txt", package="GeneGeneInteR")
dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\\t")

## In this example, genes are loosely scoured but default are much harsher
## conditions
new.snps <- snpMatrixScour(dta$snpX, genes.info = dta$genes.info,
                           min.maf = 0.2, min.eq=0.05, NA.rate = 0.2)
}
\seealso{
\code{\link{GGI}}
}

