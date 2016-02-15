#' @title Read in Lowestoft-format VPA data
#' 
#' @description Read in Lowestoft-format VPA data
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
#' @description Read Lowestoft-format survey data
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