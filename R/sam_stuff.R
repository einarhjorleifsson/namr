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

#' @title Read in Lowestoft-format VPA data
#' 
#' @export
#' 
#' @param filename The name of the file to read in
#' @param val.name XXX
#' @param format The format of the output, available is "List","Wide","Long"

read_lowestoft <- function(filename, val.name, format="List")
{
  y <- scan(filename, skip = 2, nlines = 1, quiet = TRUE)
  a <- scan(filename, skip = 3, nlines = 1, quiet = TRUE)
  tab <- read.delim(filename, header = FALSE, sep = "", skip = 5)[1:(y[2] - y[1] + 1),]
  names(tab) <- c(a[1]:a[2])
  rownames(tab) <- c(y[1]:y[2])
  
  Type <- scan(filename,skip=4,nlines=1,quiet=TRUE)
  if(Type == 2) {
   tab[2:nrow(tab),] <- tab[1,] 
  }
  
  if(format == "List") return(list(y = y, a = a, tab = tab))
  if(format == "Wide") return(tab)
  
  tab$year <- as.integer(rownames(tab))
  tab <- reshape2::melt(tab,id.vars="year",factorsAsStrings = FALSE)
  names(tab) <- c("year","age",val.name)
  tab$age <- as.integer(as.character(tab$age))
  tab <- tab[!is.na(tab$age),]
  tab[,3] <- as.numeric(tab[,3])
  return(tab)
}

#' @title read_ibya_lowestoft
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param path XXX
#' @param Scale XXX
read_ibya_lowestoft <- function(path,Scale=1) {
  oc <-  read_lowestoft(paste(path,"cn.dat",sep="/"),val.name="oC",format = "Long")
  oc$oC <- oc$oC/Scale
  cw <-  read_lowestoft(paste(path,"cw.dat",sep="/"),val.name="cW",format = "Long")
  sw <-  read_lowestoft(paste(path,"sw.dat",sep="/"),val.name="sW",format = "Long")
  mat <- read_lowestoft(paste(path,"mo.dat",sep="/"),val.name="mat",format = "Long")
  nat <- read_lowestoft(paste(path,"nm.dat",sep="/"),val.name="m",format = "Long")
  #pf <-  read_lowestoft(paste(path,"pf.dat",sep="/"),val.name="pf",format = "Long")
  #pm <-  read_lowestoft(paste(path,"pm.dat",sep="/"),val.name="pm",format = "Long")
  
  res <- plyr::join(oc,sw,by=c("year","age"))
  res <- plyr::join(res,cw,by=c("year","age"))
  res <- plyr::join(res,mat,by=c("year","age"))
  res <- plyr::join(res,nat,by=c("year","age"))
  #res <- plyr::join(res,pf,by=c("year","age"))
  #res <- plyr::join(res,pm,by=c("year","age"))
  res$oC[res$oC == -1] <- NA
  
  return(res)
}

#' @title Read Lowestoft-format survey data
#' 
#' @export
#' 
#' @param filename The name of the file to read in
#' @param format XXX

read_lowestoft_survey <- function(filename,format="long")
{
  n <- scan(filename, skip = 1, nlines = 1, quiet = TRUE) - 100
  idx <- vector("list", length = n)
  start <- 3
  for (k in 1:n)
  {
    idx[[k]]$name <- paste(scan(filename, skip = start - 1, nlines = 1, 
                                what = character(0), quiet = TRUE), collapse = " ")
    temp <- scan(filename, skip = start, nlines = 1, quiet = TRUE)
    idx[[k]]$y1 <- temp[1]
    idx[[k]]$y2 <- temp[2]
    idx[[k]]$ny <- temp[2] - temp[1] + 1
    temp <- scan(filename, skip = start + 1, nlines = 1, quiet = TRUE)
    idx[[k]]$rho <- 0.5 * (temp[4] + temp[3])
    temp <- scan(filename, skip = start + 2, nlines = 1, quiet = TRUE)
    idx[[k]]$a1 <- temp[1]
    idx[[k]]$a2 <- temp[2]
    idx[[k]]$na <- temp[2] - temp[1] + 1
    idx[[k]]$tab <- read.table(filename, skip = start + 3, nrows = idx[[k]]$ny)
    temp <- idx[[k]]$tab[,2:(idx[[k]]$na + 1)] 
    effort <- idx[[k]]$tab[,1] 
    idx[[k]]$tab <- data.frame(temp / effort)
    names(idx[[k]]$tab) <- idx[[k]]$a1:idx[[k]]$a2
    rownames(idx[[k]]$tab) <- idx[[k]]$y1:idx[[k]]$y2
    start <- start + 4 + idx[[k]]$ny
  }
  if(format=="list") return(list(n = n, idx = idx))
  
  d <- data.frame()
  
  # dummy for library check
  oU <- NULL
  
  for (i in 1:length(idx)) {
    x <- idx[[i]]$tab
    x$year <- as.integer(rownames(x))
    x <- reshape2::melt(x,"year",variable.name="age",value.name  = "oU")
    x$index <- idx[[i]]$name
    if(i == 1) {
      d <- x
    } else {
      d <- rbind(d,x)
    }
  }
  
  d$year <- as.integer(as.character(d$year))
  d$age <- as.integer(as.character(d$age))
  
  if(format == "long") {
    d <- plyr::ddply(d,c("age","index"),transform,
               oU.std=oU/mean(oU,na.rm=T),
               oU.std.log=log(oU/mean(oU,na.rm=T)))
    d$yc <- d$year - d$age
    return(d)
  }
  
  
  d <- reshape2::dcast(d,year + age ~ index,value.var = "oU")
  
  if(format == "wide") return(d)
  
  message(paste('specified',format,'not recognized. Use either "list", "long" or "wide"'))
  
  return(NULL)
}
