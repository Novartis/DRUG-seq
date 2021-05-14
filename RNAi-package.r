#' RNAi
#'
#' RNAi Package to do stuff
#'
#' @name RNAi
#'
#'

require(dplyr)



# mapper for hypergeometric distribution
# with lower.tail, return value is P[X > x], therefore x-1 as first argument
rsa.phyper.map <- function(row) {
  phyper(row[1]-1, row[2], row[3], row[4],
         lower.tail=F, log.p=T)
}

# b<-590 to 591 reaches the 1E-2180  lowest p-value limit 
# a
# phyper(b, a-b, 30203, 1090, lower.tail=FALSE, log.p = TRUE)

# calculate p-values for each siRNA of a gene
# by default, return only data for the siRNA with the smallest p-value
# with full=T, all data is returned
rsa.prob <- function(x, N, full) {
  nr <- nrow(x)
  x <- x[order(x$rank),]
  a <- seq_len(nr)
  b <- x$rank
  c <- N-b
  d <- nr
  counts <- cbind(a,b,c,d)
  p <- apply(counts, 1, rsa.phyper.map)/log(10)
  
  if (full) {
    idx <- seq_len(nr)
  } else {
    idx <- which.min(p)
  }
  
  return(cbind(x, counts, p)[idx,])
}

#' rsa
#'
#'
#' @author Florian Nigsch, Gregory McAllister
#' @name rsa
#' @param rsadat data frame containing a minimum of gene id and value
#' @param gene.column column name in rsadat specifying the gene identifier
#' @param vale.column column name in rsadat specifying the value
#' @param activator boolean specifying whether the screen is an activator or inhibitor.  Activator will ranking such that numerically smaller values get lower ranks (ie are better)
#' @param full boolean specifying whether to return all objects and their associated p-values.  This will also return counts from the hypergeometric
#' @param subset boolean (T = limit data frame returned to gene identifier and LogP)
#' @return data frame
#' @references A probability-based approach for the analysis of large-scale RNAi screens
#'   Nature Methods - 4, 847 - 849 (2007), doi:10.1038/nmeth1089
#'   Paper: \url{http://www.nature.com/nmeth/journal/v4/n10/full/nmeth1089.html}
#'   SuppInfo: \url{http://www.nature.com/nmeth/journal/v4/n10/extref/nmeth1089-S1.pdf}
#'
#' @seealso \code{rsa.gnf}
#' 
#' @export



rsa <- function(rsadat, activator=F, full=F, gene.column="gene", value.column="value", subset=F, drugseq=F,drugseq.ties='average') {
  
  switch <- ifelse(activator, -1, 1)
  
  if(drugseq) {
    ties.met<-drugseq.ties #default = 'average' which performs well upon manual review.
  } else {
    ties.met<-"min"
  }
  
  rsadat$rank <- rank(rsadat[,value.column] * switch, ties.method=ties.met)
  
  res <- by(rsadat, rsadat[,gene.column], rsa.prob, nrow(rsadat), full)
  
  res <- do.call("rbind", res)
  
  if(full)
    return(res)
  
  if (subset){
    tmp <- res[,c(gene.column,"p")]
    colnames(tmp)[2] <- "LogP"
    return(tmp)
  }else{
    tmp <- res[,!colnames(res) %in% c("rank","a","b","c","d")]
    colnames(tmp)[which(colnames(tmp)=="p")] <- "LogP"
    return(tmp)
  }
  
}


q1q3 <- function(dat){
  result <- quantile(dat$my_value,probs=c(0.25,0.75),na.rm=T)
  return (c(result[1],result[2]))
}

q0q1q3q4 <-function(dat){
  result <- quantile(dat$my_value,probs=c(0,0.25,0.75,1),na.rm=T)
  return (c(result[1],result[2],result[3],result[4]))
}

logq0q1q2q3q4<-function(dat){
  result <- quantile(dat$my_logvalue,probs=c(0,0.25,0.5,0.75,1),na.rm=T)
  return (c(result[1],result[2],result[3],result[4],result[5]))
}


#' score.screen.ZY.RQz
#'
#' Generate gene level scores for a single RNAi screen.  Runs RSA in both directions and generates Q1/Q3 per gene.
#'
#' @name score.screen
#' @param x data frame containing at least a gene column and score column
#' @param gene.column column name in rsadat specifying the gene identifier
#' @param vale.column column name in rsadat specifying the value
#' @param both.directions boolean indicating whether to calculate "RSA Up" as well
#'
#' @export

score.screen.ZY.RQz <- function(data,gene.column="gene",value.column="value",log.column="Log2FoldChange",drugseq=F,drugseq.ties='average') {
  if(drugseq){print('RSA setup for DRUG-seq: ties.methods is set to average in place of min in CRISPR/siRNA screens')}
  x.rsa.d <- rsa(data,gene.column=gene.column,value.column=value.column,subset=T,drugseq=drugseq,drugseq.ties=drugseq.ties)
  x.rsa.d <- unique(x.rsa.d[,c(gene.column,"LogP")])
  colnames(x.rsa.d)[2] <- "logP_RSA_Down"
  
  x.rsa.u <- rsa(data,gene.column=gene.column,value.column=value.column,subset=T,activator=T,drugseq=drugseq,drugseq.ties=drugseq.ties)
  x.rsa.u <- unique(x.rsa.u[,c(gene.column,"LogP")])
  colnames(x.rsa.u)[2] <- "logP_RSA_Up"
  
  #y <- ddply(data,gene.column,function(x) { quantile(x[,value.column],probs=c(0.25,0.75))})
  # remove plyr dependency - there are probably more elegant code for this thou...
  data$my_gene <- data[,gene.column]
  data$my_value <- data[,value.column]
  data$my_logvalue<- data[,log.column]
  data2<- select(data,my_gene,my_logvalue)
  data <- select(data,my_gene,my_value)
  #updated from q1q3 to q0q1q3q4 on 2017-11-20
  #y <- by(data,data$my_gene,q1q3)
  y <- by(data,data$my_gene,q0q1q3q4)
  y2<- by(data2,data2$my_gene,logq0q1q2q3q4)
  
  y <- as.data.frame(do.call(rbind,y))
  y$gene<-rownames(y)
  
  y2 <- as.data.frame(do.call(rbind,y2))
  y2$gene<-rownames(y2)
  
  #y <- y[c(5,1,2)]
  #colnames(y) <- c(gene.column,"rz_score_Q1","rz_score_Q3")
  #updated from above on 2017-11-20
  y <- y[c(5,1,2,3,4)]
  colnames(y) <- c(gene.column,"rz_score_min","rz_score_Q1","rz_score_Q3","rz_score_max")
  
  y2 <- y2[c(6,1,2,3,4,5)]
  colnames(y2) <- c(gene.column,"Min log2FoldChange","Q1 log2FoldChange","Median log2FoldChange","Q3 log2FoldChange","Max log2FoldChange")
  
  all <- merge(x.rsa.d,x.rsa.u,by=gene.column,all=T)
  
  
  zN<-by(data,data$my_gene,tally)
  zN <- as.data.frame(do.call(rbind,zN))
  zN$gene<-rownames(zN)
  zN <- zN[c(2,1)]
  colnames(zN)<-c(gene.column,"N")
  
  data_temp_D<-subset(data,data$my_value<=(-2.0))
  data_temp_D3<-subset(data,data$my_value<=(-3.0))
  data_temp_D4<-subset(data,data$my_value<=(-4.0))
  
  if(dim(data_temp_D)[1]==0) {
  } else {
    zD<-by(data_temp_D,data_temp_D$my_gene,tally)
    zD <- as.data.frame(do.call(rbind,zD))
    zD$gene<-rownames(zD)
    zD <- zD[c(2,1)]
    colnames(zD)<-c(gene.column,"z2_count_Down")
  }
  
  if(dim(data_temp_D3)[1]==0) {
  } else {
    zD3<-by(data_temp_D3,data_temp_D3$my_gene,tally)
    zD3 <- as.data.frame(do.call(rbind,zD3))
    zD3$gene<-rownames(zD3)
    zD3 <- zD3[c(2,1)]
    colnames(zD3)<-c(gene.column,"z3_count_Down")
  }
  
  if(dim(data_temp_D4)[1]==0) {
  } else {
    zD4<-by(data_temp_D4,data_temp_D4$my_gene,tally)
    zD4 <- as.data.frame(do.call(rbind,zD4))
    zD4$gene<-rownames(zD4)
    zD4 <- zD4[c(2,1)]
    colnames(zD4)<-c(gene.column,"z4_count_Down")
  }
  
  data_temp_U<-subset(data,data$my_value>=(2.0))
  data_temp_U3<-subset(data,data$my_value>=(3.0))
  data_temp_U4<-subset(data,data$my_value>=(4.0))
  
  if(dim(data_temp_U)[1]==0) {
  } else {
    zU<-by(data_temp_U,data_temp_U$my_gene,tally)
    zU <- as.data.frame(do.call(rbind,zU))
    zU$gene<-rownames(zU)
    zU <- zU[c(2,1)]
    colnames(zU)<-c(gene.column,"z2_count_Up")
  }
  
  if(dim(data_temp_U3)[1]==0) {
  } else {
    zU3<-by(data_temp_U3,data_temp_U3$my_gene,tally)
    zU3 <- as.data.frame(do.call(rbind,zU3))
    zU3$gene<-rownames(zU3)
    zU3 <- zU3[c(2,1)]
    colnames(zU3)<-c(gene.column,"z3_count_Up")
  }
  
  if(dim(data_temp_U4)[1]==0) {
  } else {
    zU4<-by(data_temp_U4,data_temp_U4$my_gene,tally)
    zU4 <- as.data.frame(do.call(rbind,zU4))
    zU4$gene<-rownames(zU4)
    zU4 <- zU4[c(2,1)]
    colnames(zU4)<-c(gene.column,"z4_count_Up")
  }
  
  all <- merge(all,y,by=gene.column,all=T)
  all <- merge(all,y2,by=gene.column,all=T)
  
  all <- merge(all,zN,by=gene.column,all=T)
  
  if(dim(data_temp_D)[1]==0) {
    all$z2_count_Down<-(0)
  } else {
    all <- merge(all,zD,by=gene.column,all=T)
  }
  
  if(dim(data_temp_U)[1]==0) {
    all$z2_count_Up<-(0)
  } else {
    all <- merge(all,zU,by=gene.column,all=T)
  }
  
  if(dim(data_temp_D3)[1]==0) {
    all$z3_count_Down<-(0)
  } else {
    all <- merge(all,zD3,by=gene.column,all=T)
  }
  
  if(dim(data_temp_U3)[1]==0) {
    all$z3_count_Up<-(0)
  } else {
    all <- merge(all,zU3,by=gene.column,all=T)
  }
  
  if(dim(data_temp_D4)[1]==0) {
    all$z4_count_Down<-(0)
  } else {
    all <- merge(all,zD4,by=gene.column,all=T)
  }
  
  if(dim(data_temp_U4)[1]==0) {
    all$z4_count_Up<-(0)
  } else {
    all <- merge(all,zU4,by=gene.column,all=T)
  }
  
  all[is.na(all$z2_count_Down),'z2_count_Down']<-0
  all[is.na(all$z2_count_Up),'z2_count_Up']<-0
  all[is.na(all$z3_count_Down),'z3_count_Down']<-0
  all[is.na(all$z3_count_Up),'z3_count_Up']<-0  
  all[is.na(all$z4_count_Down),'z4_count_Down']<-0
  all[is.na(all$z4_count_Up),'z4_count_Up']<-0
  
  invisible(all)
}