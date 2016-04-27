#'Imports SNPs information from pedfile, PLINK, VCF (4.0) file, or genotypes
#'imputed by IMPUTE2.
#'
#'\code{ImportFile} generates a
#'\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} object on the
#'basis of diallelic object contained in a file, and creates also a data frame
#'containing information about position of the SNP on the genome.
#'
#'This function uses the full path to the file to import and read the file with
#'\code{\link[snpStats]{read.pedfile}}, \code{\link[snpStats]{read.plink}},
#'\code{\link[GGtools]{vcf2sm}}, or \code{\link[snpStats]{read.impute}},
#'depending on the file extension (pedfile, plink, vcf or impute2). For a
#'pedfile, \code{ImportFile} also reads the ".info" file associated (this must
#'be in the same directory as the ".ped". Similarly, for a PLINK, there must be
#'3 files with extensions ".bed" (passed as argument \code{file}), ".bim" and
#'".fam". A VCF file must be with the associated ".tbi" file.
#'
#'If the file is a vcf file, you must precise two arguments : \itemize{\item
#'\code{gr} instance of \code{\link[GenomicRanges]{GRanges}.} \item
#'\code{nmetacol} numeric defining the number of columns used in each record as
#'locus-level metadata.}
#'
#'Pos argument is optional. If it's not given, the data frame with position
#'information is filled with NAs, except for SNP names which are imported from
#'the column names of the
#'\href{http://bioconductor.org/packages/snpStats/}{SnpMatrix}, and eventually
#'positions and chromosomes if the file imported is a pedfile or a PLINK. Else,
#'the pos argument can be either the path to a csv file, a character vector with
#'elements of the form chr:position, or a numeric vector with only the
#'positions. Additionnaly, SNP names can be precised as names of the vector. If
#'you choose the csv file path, be sure that the columns are named like
#'following : Chromosome, Genenames, SNPnames, Position. Other writings are
#'understood, refer to the function itself to know more.
#'
#'@param file String containing the path of the file to import.
#'@param pos (optional) Path to a csv file, character vector or numeric vector
#'  containing informations for each SNP about chromosome, gene names, snp names
#'  and position.
#'@param pos.sep (optional) String to be passed to \code{read.csv} function as
#'  \code{sep} argument. Default is tab (\\t).
#'@param ... Additionnal arguments to be passed to the reading file function.
#'@return The output is a list of two objects : \code{snpX}, a
#'  \href{http://bioconductor.org/packages/snpStats/}{SnpMatrix} and \code{genes.info}
#'  a data frame with 4 columns, and one row per SNP. The columns are
#'  Chromosome, Genenames, SNPnames and position.
#'@export
#'
#'@importFrom IRanges IRanges
#'
#'@importFrom GenomicRanges GRanges
#'
#' @examples
#' ## Pedfile from this package.
#' ped <- system.file("extdata/example.ped", package="GGItest")
#' info <- system.file("extdata/example.info", package="GGItest")
#'
#' ## Information about position of the snps
#' posi <- system.file("extdata/example.txt", package="GGItest")
#'
#' ## Importation
#' dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\t")
#'
#' #########
#' ## VCF file from GGtools package.
#' vref <- system.file("vcf/CEU.exon.2010_09.genotypes.vcf.gz", package="GGtools")
#' irange <- IRanges::IRanges(10e6,20e6)
#' gg = GenomicRanges::GRanges(seqnames="1", ranges=irange)
#' dta <- ImportFile(file=vref, gr=gg, nmetacol=9L)

ImportFile <- function (file, pos, pos.sep="\t", ...) {

  if(class(file)!="character"){
    stop("Please, enter a valid path file for genotype data.")
  }

  nc <- nchar(file)
  ext <- substr(file, nc - 3, nc)
  extgz <- substr(file, nc - 6, nc)
  res <- list()
  imp <- NULL

  if (ext==".ped") {
    imp <- snpStats::read.pedfile(file = file, ...)

  } else if (extgz==".ped.gz") {
    imp <- snpStats::read.pedfile(file = file, ...)

  } else if (ext == ".bed") {
    imp <- snpStats::read.plink(bed = file, ...)

  } else if (ext %in% c(".vcf", "f.gz")) {
    imp <- GGtools::vcf2sm(tbxfi = Rsamtools::TabixFile(file), ...)

  } else if (ext == ".impute2") {
    imp <- snpStats::read.impute(file, ...)

  } else {stop("Please enter a valid pedfile, plink, vcf, or impute2 file.")}

  if(class(imp)=="list"){
    res[["snpX"]] <- imp$genotypes

  } else if(class(imp)=="SnpMatrix"){
    res[["snpX"]] <- imp
  }

  if(!missing(pos)){

    if(length(pos)==1){
      if(class(pos)=="character"){
        infos <- read.csv(pos, sep=pos.sep, header=T)

        chr <- infos[,names(infos)%in%c("Chromosome","Chr","chromosome","chr")]
        gene <- infos[,names(infos)%in%c("Gene","gene","genenames","Genenames","Gene.names","gene.names")]
        snp <- infos[,names(infos)%in%c("SNP","Snp","snp","SNPnames","Snpnames","snpnames","SNP.names","Snp.names","snp.names")]
        posi <- infos[,names(infos)%in%c("Position","position","pos","Pos")]

        if(length(chr)==0){chr<-rep(NA,nrow(infos));warning("Chromosome column was not found.")}
        else if(length(gene)==0){gene<-rep(NA,nrow(infos));warning("Gene names column was not found.")}
        else if(length(snp)==0){snp<-colnames(res[["snpX"]])}
        else if(length(posi)==0){posi<-rep(NA,nrow(infos));warning("Position column was not found.")}

        res[["genes.info"]] <- data.frame(Chromosome=chr,
                                    Genenames=gene,
                                    SNPnames=snp,
                                    Position=posi)
      } else {
        stop("Pos argument needs to be either a numeric vector, a character vector or a path file.")
      }
    } else if(length(pos)==ncol(res[["snpX"]])){
      if(class(pos)%in%c("numeric","integer")){
        snp <- names(pos)
        if(is.null(snp)){snp <- colnames(res[["snpX"]])}
        res[["genes.info"]] <- data.frame(Chromosome=rep(NA,length(pos)),
                                    Genenames=rep(NA,length(pos)),
                                    SNPnames=snp,
                                    Position=pos)

      } else if(class(pos)=="character"){
        snp <- names(pos)
        if(is.null(snp)){snp <- colnames(res[["snpX"]])}

        liste <- data.table::tstrsplit(pos, split=":", fixed=TRUE)
        res[["genes.info"]] <- data.frame(Chromosome=liste[[1]],
                                    Genenames=rep(NA,length(pos)),
                                    SNPnames=snp,
                                    Position=liste[[2]])
      } else {
        stop("Pos argument needs to be either a numeric vector, a character vector or a path file.")

      }
    } else {
      stop("The number of SNPs must be the same in genotype data and position information.")
    }
  } else {
    if(class(imp)=="list"){
      chr <- imp$map[,names(imp$map)%in%c("Chromosome","Chr","chromosome","chr")]
      gene <- imp$map[,names(imp$map)%in%c("Gene","gene","genenames","Genenames","Gene.names","gene.names")]
      snp <- imp$map[,names(imp$map)%in%c("SNP","Snp","snp","SNPnames","Snpnames","snpnames","SNP.names","Snp.names","snp.names","snp.name","SNP.name","Snp.name","SNPname","Snpname","snpname")]
      posi <- imp$map[,names(imp$map)%in%c("Position","position","pos","Pos","V2")]

      if(class(posi)=="character"){
        liste <- data.table::tstrsplit(posi, split=":", fixed=TRUE)
        chr <- liste[[1]]
        posi <- liste[[2]]
      }

      if(length(chr)==0){chr<-rep(NA,nrow(imp$map));warning("Chromosome column was not found.")}
      if(length(gene)==0){gene<-rep(NA,nrow(imp$map));warning("Gene names column was not found.")}
      if(length(snp)==0){snp<-colnames(res[["snpX"]])}
      if(length(posi)==0){posi<-rep(NA,nrow(imp$map));warning("Position column was not found.")}

      res[["genes.info"]] <- data.frame(Chromosome=chr,
                                  Genenames=gene,
                                  SNPnames=snp,
                                  Position=posi)
    } else {
      chr<-rep(NA,ncol(res[["snpX"]]))
      gene<-rep(NA,ncol(res[["snpX"]]))
      snp<-colnames(res[["snpX"]])
      posi<-rep(NA,ncol(res[["snpX"]]))

      res[["genes.info"]] <- data.frame(Chromosome=chr,
                                  Genenames=gene,
                                  SNPnames=snp,
                                  Position=posi)
    }
  }

  if(any(colnames(res[["snpX"]])!=res[["genes.info"]][,"SNPnames"])){warning("Be careful, the SNP names don't match between snpMatrix and info dataframe.")}

  return(res)
}
