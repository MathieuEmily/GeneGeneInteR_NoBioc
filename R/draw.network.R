draw.network <- function(GGI,genes=1:ncol(GGI),threshold=0.05,plot.nointer=TRUE){

  if(length(genes)<2 || length(genes)>ncol(GGI)){
    stop("Number of genes selected not valid.")
  } else if(!class(GGI)%in%c("data.frame","matrix")){
    stop("GGI must be a data.frame.")
  } else if(ncol(GGI)!=nrow(GGI)){
    stop("GGI must be a sqared matrix, containing the pValues for each interaction between genes.")
  } else if(!class(threshold)%in%c("numeric","integer")){
    stop("Threshold must be a numeric.")
  } else if(threshold>1 || threshold<0){
    stop("Threshold must be comprised in [0,1].")
  } else if(class(plot.nointer)!="logical"){
    stop("plot.inter must be a boolean.")
  }

  if(class(genes)=="character"&&any(!genes%in%colnames(GGI))){
    stop("Genes and GGI don't match. Please select genes that are named in GGI.")
  }

  GGI <- GGI[genes,genes]
  dim <- ncol(GGI)
  pVal.raw <- GGI[lower.tri(GGI)]
  if(any(is.na(pVal.raw))){
    warning("NAs found in GGI, considered as not significative.")
    pVal.raw[is.na(pVal.raw)]<-1
  }

  from.raw <- c()
  to.raw <- c()

  for (i in 1:(dim-1)){
    from.raw <- c(from.raw, rep(colnames(GGI)[i], dim-i))
    to.raw <- c(to.raw, rownames(GGI)[(i+1):dim])
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

}
