#'Draw a gene interactions network.
#'
#'\code{\link{draw.network}} generates graph representing the significative
#'interactions between genes of a GGI study, using the output of
#'\code{\link{GGI}}.
#'
#'The function takes as input the output of \code{GGI}, a matrix of pvalues
#'representing the level of interaction between genes. \code{genes} can be used
#'to select the genes you want to plot on the interaction network. It is either
#'a character vector giving genes' name or a numeric vector giving their column
#'number.
#'
#'@param gene.interactions Data.frame of pValues containing the results of a GGI
#'  analysis. The data.frame must be a squared matrix.
#'@param genes Numeric vector allowing a selection of the genes that will be
#'  included in the relations. Default is set to all genes.
#'@param threshold Numeric comprised in [0,1], defining the pValues that will be
#'  considered as respresentative of an interaction. Default is set to 0.05.
#'@param plot.nointer A boolean. Set TRUE if the genes with no interaction
#'  should be plotted.
#'@return The output is a list of three objects : \code{from} and \code{to},
#'  that are character vectors giving the genes that interact, and \code{pVal},
#'  a numeric vector giving the pvalues for each interaction. A warning message
#'  is displayed if the vector \code{from} is empty, which means no interaction
#'  has been found, according to the threshold given.
#'@export
#'
#'@seealso \code{\link{draw.network}}, \code{\link{GGI}}
#'
#' @examples
#'
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
#' pVal <- GGI(Y=Y, snpX=dta$snpX, genes.info=dta$genes.info)
#'
#' ## By default, all genes are selected, threshold is 0.05.
#' draw.network(pVal)
#'
#' ## By default, genes with no interaction are plotted, use plot.nointer to change this option.
#' draw.network(pVal, threshold=0.005)
#' draw.network(pVal, threshold=0.005, plot.nointer=F)

draw.network <- function(gene.interactions,genes=1:ncol(gene.interactions),threshold=0.05,plot.nointer=T){

  if(length(genes)<2 || length(genes)>ncol(gene.interactions)){
    stop("Number of genes selected not valid.")
  } else if(!class(gene.interactions)%in%c("data.frame","matrix")){
    stop("Gene.interactions must be a data.frame.")
  } else if(ncol(gene.interactions)!=nrow(gene.interactions)){
    stop("Gene.interactions must be a sqared matrix, containing the pValues for each interaction between genes.")
  } else if(!class(threshold)%in%c("numeric","integer")){
    stop("Threshold must be a numeric.")
  } else if(threshold>1 || threshold<0){
    stop("Threshold must be comprised in [0,1].")
  } else if(class(plot.nointer)!="logical"){
    stop("plot.inter must be a boolean.")
  }

  if(class(genes)=="character"&&any(!genes%in%colnames(gene.interactions))){
    stop("Genes and gene.interactions don't match. Please select genes that are named in gene.interactions.")
  }

  gene.interactions <- gene.interactions[genes,genes]
  dim <- ncol(gene.interactions)
  pVal.raw <- gene.interactions[lower.tri(gene.interactions)]
  if(any(is.na(pVal.raw))){
    warning("NAs found in gene.interactions, considered as not significative.")
    pVal.raw[is.na(pVal.raw)]<-1
  }

  from.raw <- c()
  to.raw <- c()

  for (i in 1:(dim-1)){
    from.raw <- c(from.raw, rep(colnames(gene.interactions)[i], dim-i))
    to.raw <- c(to.raw, rownames(gene.interactions)[(i+1):dim])
  }



  from <- from.raw[pVal.raw<threshold]
  to <- to.raw[pVal.raw<threshold]
  pVal <- pVal.raw[pVal.raw<threshold]

  if(plot.nointer){
    actors <- data.frame(name=levels(as.factor(unique(c(from.raw,to.raw)))))
  } else {
    actors <- data.frame(name=levels(as.factor(unique(c(from,to)))))
  }

  if(length(from)==0){warning("No interactions has been found between the genes selected.")}

  relations <- data.frame(from=from,
                          to=to,
                          pVal=pVal)

  g <- igraph::graph_from_data_frame(relations, directed=F, vertices=actors)
  plot(g, vertex.size=10)

  return(NULL)
}
