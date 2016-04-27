#'Selects SNPs in snpXatrix and associated info data.frame, outputs of
#'ImportFile.
#'
#'\code{select.snps} allows the user to select snps from an object output of
#'\code{\link{ImportFile}}. It generates the same object, with the columns of
#'the \href{http://bioconductor.org/packages/snpStats/}{snpXatrix} and the rows
#'of the data.frame corresponding to the selected snps.
#'
#'For a good working of the function, let the column names of the
#'\code{genes.info} data.frame as they were when outputted from
#'\code{\link{ImportFile}}.
#'
#'\code{select} argument can be either : \itemize{\item a numeric vector with
#'only the column number in the
#'\href{http://bioconductor.org/packages/snpStats/}{snpXatrix} (or row number
#'for \code{genes.info}) of each snp selected.\item a character vector with the
#'names of each snp selected or each gene selected.\item a character vector
#'which elements are position bounds of genes. Each element of the vector is
#'either of the form "begin:end", or "chr:begin:end" if you have to precise the
#'chromosome of the gene.}
#'
#'@param snpX \href{http://bioconductor.org/packages/snpStats/}{snpXatrix}
#'  object. Given as output of \code{\link{ImportFile}}.
#'@param genes.info Data.frame containing informations about snps. For more
#'  details, refer to \code{\link{ImportFile}} help file.
#'@param select Numeric or character vector for selecting snps in \code{snpX}
#'  and \code{genes.info}. See details for more information.
#'@return The output is a list of two objects : \code{snpX}, a
#'  \href{http://bioconductor.org/packages/snpStats/}{snpXatrix} and \code{genes.info}
#'  a data frame with 4 columns, and one row per SNP selected with
#'  \code{select}. The columns are Chromosome, Genenames, SNPnames and Position.
#'  An error message is displayed if the genes of snps selected are not found in
#'  the either \code{snpX} or \code{genes.info}.
#'@export
#'
#' @examples
#'
#' ## Dataset is included in the package
#' ped <- system.file("extdata/example.ped", package="GGItest")
#' info <- system.file("extdata/example.info", package="GGItest")
#' posi <- system.file("extdata/example.txt", package="GGItest")
#' dta <- ImportFile(file=ped, snps=info, pos=posi, pos.sep="\t")
#'
#' ## Selection of the genes DNAH9 and TXNDC5
#' selec <- select.snps(dta$snpX, dta$genes.info, c("DNAH9","TXNDC5"))
#'
#' #' ## Selection of the snps from position 101342000 to 101490000 on chromosome 15
#' selec <- select.snps(dta$snpX, dta$genes.info, c("15:101342000:101490000"))

select.snps <- function(snpX, genes.info, select){
  if(class(snpX)!="SnpMatrix"){
    stop("snpX must be a snpMatrix object.")
  } else if(class(genes.info)!="data.frame"){
    stop("genes.info must be a data.frame.")
  } else if(ncol(snpX)!=nrow(genes.info)){
    stop("snpX and genes.info must have the same number of snps.")
  }

  res <- NULL

  if(class(select)%in%c("numeric","integer")){
    if(any(c(select<1, select>ncol(snpX)))){
      stop("The selection of snps is out of bounds.")
    } else {
      res[["gentoypes"]] <- snpX[,select]
      res[["genes.info"]] <- genes.info[select,]
    }

  } else if(class(select)=="character"){
    if(all(select %in% colnames(snpX))){
      res[["snpX"]] <- snpX[,select]
      res[["genes.info"]] <- genes.info[genes.info[,"SNPnames"]%in%select,]
    } else if(all(select %in% genes.info[,"Genenames"])){
      genes <- genes.info[,"Genenames"]%in%select
      res[["snpX"]] <- snpX[,as.character(genes.info[genes,"SNPnames"])]
      res[["genes.info"]] <- genes.info[genes,]
    } else {
      liste <- strsplit(select, ":")
      snp <- c()
      if(length(liste[[1]])==2){
        for(i in 1:length(liste)){
          rows <- genes.info[as.numeric(genes.info[,"Position"])>=as.numeric(liste[[i]][1]),]
          rows <- rows[as.numeric(rows[,"Position"])<=as.numeric(liste[[i]][2]),]
          snp <- c(snp,as.character(rows[,"SNPnames"]))
        }
      } else if(length(liste[[1]])==3){
        for(i in 1:length(liste)){
          rows <- genes.info[as.numeric(genes.info[,"Chromosome"])==as.numeric(liste[[i]][1]),]
          rows <- rows[as.numeric(rows[,"Position"])>=as.numeric(liste[[i]][2]),]
          rows <- rows[as.numeric(rows[,"Position"])<=as.numeric(liste[[i]][3]),]
          snp <- c(snp,as.character(rows[,"SNPnames"]))
        }
      } else {stop("Wrong format for select argument, or genes or snps not found. Please refer to help file.")}

      if(length(snp)==0){stop("No matches were found with your selection.")}
      else{
        res[["snpX"]] <- snpX[,snp]
        res[["genes.info"]] <- genes.info[genes.info[,"SNPnames"]%in%snp,]
      }
    }

  } else {
    stop("The argument select must be a numeric or character vector.")
  }

  rownames(res[["genes.info"]]) <- 1:nrow(res[["genes.info"]])
  res$genes.info <- droplevels(res$genes.info)
  return(res)
}
