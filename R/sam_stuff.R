#' @title read_conf_sam
#' 
#' @description Reads sam configuration file and returns a list. 
#' Used by read_rby_sam, xxx.
#' This function is identical(?) to originally in scr/common.R 
#' 
#' @author Anders Nielsen
#' 
#' @export
#' 
#' @param filename A character vector containg path and name to ....

read_conf_sam <- function(filename){
  lin<-readLines(filename)
  idxNam<-grep("^[[:alpha:]]",lin)
  doone<-function(i){
    idxNam<-c(idxNam,length(lin)+1)
    if(abs(idxNam[i]-idxNam[i+1])>2){
      x<-read.table(textConnection(lin[(idxNam[i]+1):(idxNam[i+1]-1)]))
      names(x)<-NULL
      as.matrix(x)
    }
  }
  ret<-lapply(1:length(idxNam),doone)
  names(ret)<-sub(' =','',lin[idxNam])
  return(ret)
}

#' @title read_rby_sam
#' 
#' @description Read sam results by year and age
#' 
#' @export
#' 
#' @param conf_sam An list object returned by function \code{read_conf_sam}
#' @param format A character vector specifying output format. If default ("long") then function
#' returns a data.frame ....

read_rby_sam <- function(conf_sam,format="long") {
  
  x <- conf_sam
 
  if(format == "wide") {
    d <- cbind(x$year,x$R[,-2],x$tsb[,-2],x$ssb[,-2],x$fbar[,-2],exp(x$logCatch[,-2]))
    colnames(d) <- c("year","r","r1","r2","b","b1","b2","ssb","ssb1","ssb2",
                   "f","f1","f2","y","y1","y2")
    d <- as.data.frame(d)
    return(d)
  }
  
  if (format == "long") {
    d <- cbind(x$year,x$R[,-2])
    colnames(d) <- c("year","est","low","hig")
    d <- as.data.frame(d)
    d$par <- "r"
    #d <- reshape2::melt(d,c("year","par"))
    
    y <- as.data.frame(cbind(x$year,x$tsb[,-2]))
    names(y)[1] <- "year"
    y$par <- "bio"
    d <- rbind(d,y)
    #d <- rbind(d,reshape2::melt(y,c("year","par")))
    
    y <- as.data.frame(cbind(x$year,x$ssb[,-2]))
    names(y)[1] <- "year"
    y$par <- "ssb"
    d <- rbind(d,y)
    #d <- rbind(d,reshape2::melt(y,c("year","par")))
    
    y <- as.data.frame(cbind(x$year,x$fbar[,-2]))
    names(y)[1] <- "year"
    y$par <- "f"
    d <- rbind(d,y)
    #d <- rbind(d,reshape2::melt(y,c("year","par")))
    
    y <- as.data.frame(cbind(x$year,exp(x$logCatch[,-2])))
    names(y)[1] <- "year"
    y$par <- "catch"
    d <- rbind(d,y)
    #d <- rbind(d,reshape2::melt(y,c("year","par")))
    
    return(d)
  }
  
}

#' @title read_fit_sam
#' 
#' @description Reads the stock assessment summary statistics from the
#' stockassessment.org directory structure
#' 
#' This function originally in scr/common.R as read.fit
#' 
#' @export
#' 
#' @param file Path to the directory containing sam run. Normally ends with
#' ../run/sam
#' @param reduced Boolean, if FALSE (default) ....
#' 
#' @return Returns a list
#' 
#' @note TO DO: Ideally should convert the data to FLStock object and from there
#' convert the stuff to e.g. rby or rbya.
#' 
#' TO DO: Provide detail description of the list returned
#' 
read_fit_sam <- function(file, reduced=FALSE){
  
  # Function to read a basic fit
  ret <- list()
  parfile <- as.numeric(scan(paste(file,'.par', sep=''), 
                             what='', n=16, quiet=TRUE)[c(6,11,16)])
  ret$nopar <- as.integer(parfile[1])
  ret$nlogl <- parfile[2]
  ret$maxgrad <- parfile[3]
  rep<-scan(paste(file,'.rep', sep=''), quiet=TRUE)
  ret$res<-read.table(paste(file,'.res', sep=''),header=FALSE)
  ret$stateDim<-rep[1]
  ret$years<-rep[-1]
  file<-paste(file,'.cor', sep='')
  lin<-readLines(file)
  ret$npar<-length(lin)-2
  ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2])
  sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!=''])
  ret$names<-unlist(lapply(sublin,function(x)x[2]))
  ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3])))
  ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4])))
  
  ret$cor<-matrix(NA, ret$npar, ret$npar)
  for(i in 1:ret$npar){
    ret$cor[1:i,i]<-as.numeric(unlist(lapply(sublin[i],
                                             function(x)x[5:(4+i)])))
    ret$cor[i,1:i]<-ret$cor[1:i,i]
  }
  ret$cov<-ret$cor*(ret$std%o%ret$std)
  
  mslh <- function(name){
    idx<-which(ret$names==name)
    x<-cbind(ret$est[idx], ret$std[idx], ret$est[idx]-2*ret$std[idx], 
             ret$est[idx]+2*ret$std[idx])
    colnames(x)<-c('est', 'std', 'low', 'hig')
    return(x)
  }
  
  ret$ssb<-mslh('ssb')
  ret$fbar<-mslh('fbar')
  ret$tsb<-mslh('tsb')
  ret$logssb<-mslh('logssb')
  ret$logfbar<-mslh('logfbar')
  ret$logtsb<-mslh('logtsb')
  ret$logscale<-mslh('logScale')
  ret$logFpar<-mslh('logFpar')
  ret$logCatch<-mslh('logCatch')
  x<-mslh('U')
  ret$stateEst<-matrix(x[,1],ncol=ret$stateDim, byrow=TRUE)
  ret$stateStd<-matrix(x[,2],ncol=ret$stateDim, byrow=TRUE)
  ret$stateLow<-matrix(x[,3],ncol=ret$stateDim, byrow=TRUE)
  ret$stateHig<-matrix(x[,4],ncol=ret$stateDim, byrow=TRUE)
  ret$R<-cbind(exp(ret$stateEst[,1]), NA, exp(ret$stateLow[,1]), 
               exp(ret$stateHig[,1]))
  if(reduced){
    ret <- ret[which(!names(ret)%in%c('cov','cor'))]
  }
  
  file <- sub('[[:alpha:]]+\\.cor$','confclone.log',file)
  if(file.exists(file)){
    ret$keys <- read_conf_sam(file)
  }
  return(ret)
}

#' @title read_rbya_sam
#' 
#' @description read fishing mortality and stock in numbers from sam
#' 
#' @export 
#' 
#' @param x Object from read_fit_sam
#' @param Scale A value 
#' 
read_rbya_sam <- function(x,Scale=1) {
  
  minAge <- min(x$res[,3])  
  maxAge <- max(x$res[,3][x$res[,3]<98.5])  
  noN <- maxAge - minAge+1
  noFleet <- max(x$res[,2])
  
  N <- exp(x$stateEst[,1:noN])/Scale
  colnames(N) <- c(minAge:maxAge)
  rownames(N) <- x$years
  N <- reshape2::melt(N,factorsAsStrings = FALSE)
  
  mort <- exp(x$stateEst[,-c(1:noN)])[,x$keys$keyLogFsta[1,]]
  colnames(mort) <- c(minAge:maxAge)
  rownames(mort) <- x$years
  mort <- reshape2::melt(mort,factorsAsStrings = FALSE)
  
  res <- cbind(N,mort[,3])
  names(res) <- c("year","age","n","f")
  
  
  # Here we need to read in the input file by read_ibya_lowestoft
  # And join that with the results.
  return(res)
}


