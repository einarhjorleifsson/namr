#' @title read_result_hake
#' 
#' @description Reads the ADMB report (result) file from the Namibian hake assessment
#' 
#' @export
#' 
#' @param infilename A character string incuding the path (directory of the input file name.
#' For hake the name of the file is normally namhakdata.dat
#' @param outfilename A character string, including path (directory) to the report
#' file. The report file has normally the ending .rep.
#' @param runname A character specifying the name of the run

read_rbx_hake <- function(infilename,outfilename,runname) {
  
  txt <- readLines(outfilename)
  txt <- stringr::str_trim(txt)
  
  tmpfile <- tempfile()
  
  # SigmaR_input SigmaR_output MSY Bsp/K
  p1 <- grep("SigmaR_input",txt)
  x <- str_clean(txt,p1,p1+1,header=T)
  x <- reshape2::melt(x)
  names(x) <- c("V1","V2")
  
  # Alpha_first
  p1 <- grep("Alpha_first",txt)
  p2 <- p1
  g <- txt[(p1):(p2)]
  g <- stringr::str_trim(g)
  g <- unlist(strsplit(g," "))
  x2_names <- g[c(1,3,5,7,9)]
  x2_values <- as.numeric(g[c(2,4,6,8,10)])
  x2 <- data.frame(V1=x2_names,V2=x2_values)
  sr_parameters <- x2
  names(sr_parameters) <- c("variable","value")
  x <- rbind(x,x2)
  
  # 'ln(Likelihoods)
  p1 <- grep("Likelihoods",txt) + 1
  p2 <- grep("M                  ",txt)
  x2 <- str_clean(txt,p1,p2)
  
  likelihoods <- x2[1:9,]
  names(likelihoods) <- c("variable","value")
  
  refpts <- x2[c(10,11,15,16,17),]
  names(refpts) <- c("variable","value")
  refpts <- rbind(refpts,data.frame(variable="Fmsy",value=refpts$value[5]/refpts$value[4]))
  refpts$details <- c("Virgin ssb","Virgin exploitable biomass",
                      "SSB at MSY","Exploitable biomass at MSY",
                      "Maximum sustainable yield",
                      "Fishing mortality resulting in MSY")
  
  x <- rbind(x,x2)
  
  p1 <- grep("sigmaCPUE1",txt)
  p2 <- grep("Addvariance",txt)
  x2 <- str_clean(txt,p1,p2)
  x <- rbind(x,x2)
  
  p1 <- grep("sigCAA_com",txt)
  p2 <- grep("sigSeal_Index",txt)
  x2 <- str_clean(txt,p1,p2)
  x <- rbind(x,x2)
  
  p1 <- grep("Natural mortality by age",txt) + 1
  p2 <- grep("M\\(8\\)",txt)
  x2 <- str_clean(txt,p1,p2,header=F)
  x <- rbind(x,x2)
  
  p1 <- grep("CPUEseries1",txt)
  p2 <- grep("Survseries2",txt)
  x2 <- str_clean(txt,p1,p2,header=F)
  x2$V1 <- paste(x2$V1,x2$V2,sep="_")
  x2$V4 <- paste(x2$V1,x2$V4,sep="_")
  x2a <- x2[,c(1,3)]
  names(x2a) <- c("V1","V2")
  x <- rbind(x,x2a)
  x2a <- x2[,c(4,5)]
  names(x2a) <- c("V1","V2")
  x <- rbind(x,x2a)
  
  x <- x[!duplicated(x),]
  names(x) <- c("variable","value") 
  tmp_out <- x
  
  # INDICES - the first results by year
  p1 <- grep('Year obs1 pred1 obs2 pred2 obs3 pred3 obs4 pred4 obs5',txt)
  p2 <- grep('Virgin Age Structure',txt)
  g <- txt[(p1):(p2-2)]
  # amend the first line (header)
  g[1] <- paste('year','oU1','pU1','oU2','pU2','oU3','pU3','oU4','pU4','oU5','pU5',
                'oU6','pU6','oUs','pUs','oUw','pUw','oUss','xx','pUss',sep=' ')
  rby <- str_clean(g,1,length(g),sep=' ',header=T)
  # The second last colum (now labelled 'xx' is redundant). Lets get rid of it:
  rby$xx <- NULL
  # ZERO's are actually NA - get rid of them
  rby <- reshape2::melt(rby,id.vars='year')
  rby$value[rby$value == 0] <- NA
  rby <- reshape2::dcast(rby,year ~ variable,value.var='value')
  
  # STOCK IN NUMBERS
  p1 <- grep('Numbers at Age',txt) + 1
  p2 <- grep('Year Bsp Bexp F RecRes Catch Depletion TotalB',txt) - 2
  rbya <- str_clean(txt,p1,p2,header=F)
  names(rbya) <- c('year',0:8)
  rbya <- reshape2::melt(rbya,id.vars='year')
  names(rbya) <- c('year','age','n')
  rbya$age <- as.integer(as.character(rbya$age))
  
  # SOME MORE results by year
  p1 <- grep('Year Bsp Bexp F RecRes Catch Depletion TotalB',txt)
  p2 <- grep('Commercial selectivity',txt)
  tmp <- str_clean(txt,p1,p2-2,header=T)
  names(tmp) <- c('year','ssb','biofish','f','rR','y','dep','bio')
  rby <- plyr::join(rby,tmp,by="year")
  
  # COMMERCIAL SELECTION PATTERN
  p1 <- p2 + 1
  p2 <- grep("Observed catch in numbers by age",txt) - 2
  sel <- str_clean(txt,p1,p2,header=T,fill=T)
  names(sel) <- c('year',0:8)
  sel <- reshape2::melt(sel,id.vars='year')
  names(sel) <- c('year','age','sel')
  sel$age <- as.numeric(as.character(sel$age))
  rbya <- plyr::join(rbya,sel,by=c("year","age"))
  
  # HERE WILL SKIP READING SOME FILES - BECAUSE SAME INFO COMES LATER IN THE OUTPUT
  
  p1 <- grep('Replacement Yield by year',txt)
  p2 <- grep('----------Stock-recruitment curve',txt)
  x <- str_clean(txt,p1+1,p2-1,header=F)
  names(x) <- c('year','rY','bio','landings') # some guesswork
  x <- x[,c("year","rY","landings")]
  rby <- plyr::join(rby,x,by="year")
  
  # Add recruitment to rby
  rby$rec <- rbya[rbya$age == 0,"n"]
  
  # Reorder column - just for aesthetic reasons
  rby <- rby[,c("year","rec","bio","ssb","biofish","f","landings","rY","dep",
                "oU1","pU1","oU2","pU2","oU3","pU3","oU4","pU4","oU5","pU5",
                "oU6","pU6","oUs","pUs","oUw","pUw","oUss","pUss")]
  rby$landings[rby$landings == 0] <- NA
  rby$f[rby$f == 0] <- NA
  # Create a separate object for the observed vs fitted value
  x <- rby[,c(1,10:ncol(rby))]
  x.obs <- x[,c(1,2,4,6,8,10,12,14,16,18)]
  x.pre <- x[,c(1,3,5,7,9,11,13,15,17,19)]
  cn <- c("year","ICSEAF (1.3 + 1.4)","ICSEAF -1.5","GLM","winter Spanish",
          "summer Spanish","7vessel CPUE","summer","winter","seal scat")  	
  names(x.obs) <- names(x.pre) <- cn
  x.obs <- reshape2::melt(x.obs,"year")
  names(x.obs) <- c("year","fleet","obs")
  x.pre <- reshape2::melt(x.pre,"year")
  names(x.pre) <- c("year","fleet","pre")
  rby.fit <- plyr::join(x.obs,x.pre,by=c("year","fleet"))
  
  # STOCK RECRUITMENT STUFF
  p1  <- grep('----------Stock-recruitment curve',txt)
  p2 <- grep('----------CommercialCAA',txt)
  xy_ssb.r <- str_clean(txt,p1+1,p2-2,header=F)
  xy_ssb.r <- as.data.frame(t(xy_ssb.r))
  names(xy_ssb.r) <- c('ssb','r')
  
  #### Catch at age
  p1 <- grep('----------CommercialCAA',txt)
  p2 <- grep('----------Survey1',txt)[2]
  x <- str_clean(txt,p1+1,p2-1)
  names(x) <- c('year','age','oC','pC','rC')
  rbya <- plyr::join(rbya,x,by=c("year","age"))
  
  # Summer survey at age
  p1 <- grep('----------Survey1',txt)[2]
  p2 <- grep('----------Survey2',txt)[2]
  x <- str_clean(txt,p1+1,p2-1,header=F)
  names(x) <- c('year','age','oU1','pU1','rU1')
  rbya <- plyr::join(rbya,x,by=c("year","age"))
  
  #### Winter survey residuals
  p1 <- grep('----------Survey2',txt)[2]
  p2 <- grep('biomass of fish by age',txt)
  x <- str_clean(txt,p1+1,p2-2,header=F)
  names(x) <- c('year','age','oU2','pU2','rU2')
  rbya <- plyr::join(rbya,x,by=c("year","age"))
  
  #### Biomass at age
  p1 <- grep('biomass of fish by age',txt)
  p2 <- grep('biomass of fish that died due to natural mortality by age for J-P Roux',txt)
  x <- str_clean(txt,p1+1,p2-2,header=F)
  g <- txt[(p1+1):(p2-2)]
  writeLines(g,tmpfile)
  x <- read.table(tmpfile,header=FALSE)
  names(x) <- c('year',0:8)
  x <- reshape2::melt(x,id.vars='year')
  names(x) <- c('year','age','b')
  x$age <- as.integer(as.character(x$age))
  rbya <- plyr::join(rbya,x,by=c("year","age"))
  
  if(!missing(runname)) {
    rby$name <- runname
    rbya$name <- runname
    xy_ssb.r$name <- runname
    tmp_out$Name <- runname
  }
  
  ## Extract some results by age data
  # Natural mortality
  i <- tmp_out$variable %in% c("M(0)","M(1)","M(2)","M(3)","M(4)","M(5)","M(6)","M(7)","M(8)")
  rba <- tmp_out[i,]
  rba$age <- rep(0:8)
  rba <- rba[,c("age","value")]
  names(rba) <- c("age","m")
  
  rbya$m <- rep(rba$m,(max(rbya$year)-min(rbya$year)+1))
  
  # quick fix here
  x <- read_input_hake(infilename)
  names(x) <- c("year","age","sW","cW","mat")
  rbya <- plyr::join(rbya,x,by=c("year","age"))
  i <- x$year == min(x$year)
  x <- x[i,]
  rba$sW <- x$sW
  rba$cW <- x$cW
  rba$mat <- x$mat
  rbya$yc <- rbya$year - rbya$age
  
  rownames(rba) <- NULL
  
  rbx <- list(rby=rby,rba=rba,rby.fit=rby.fit,rbya=rbya,xy_ssb.r=xy_ssb.r,
              likelihoods=likelihoods,refpts=refpts,tmp_out=tmp_out)
  
  return(rbx)
}


#' @title read_input_hake
#' 
#' @description Reads the ADMB input data file from the Namibian hake assessment
#' 
#' @export
#' 
#' @param filename A character string, including path (directory) to the report
#' file. The report file has normally the ending .rep.

read_input_hake <- function(filename) {
  
  tmpfile <- tempfile()
  
  txt <- readLines(filename)
  txt <- stringr::str_trim(txt)
  
  # Weight-at-age  
  # Start-Year																													
  #	0	1	2	3	4	5	6	7	8	
  p1 <- grep('Weight-at-age',txt)
  p2 <- grep('Mid-Year',txt)
  g <- txt[(p1+3):(p2-1)]
  g <- stringr::str_trim(g)
  g <- stringr::str_replace_all(g,"\\s+","\t")
  #g <- as.data.frame(g)
  writeLines(g,tmpfile)
  w1 <- as.data.frame(read.table(tmpfile,header=FALSE))
  
  txt <- readLines(filename)
  txt <- stringr::str_trim(txt)
  
  # Weight-at-age  
  # Start-Year  																												
  #	0	1	2	3	4	5	6	7	8	
  p1 <- grep('Mid-Year',txt)
  p2 <- grep('Catches',txt)
  g <- txt[(p1+2):(p2-1)]
  g <- stringr::str_trim(g)
  g <- stringr::str_replace_all(g,"\\s+","\t")
  #g <- as.data.frame(g)
  writeLines(g,tmpfile)
  w2 <- as.data.frame(read.table(tmpfile,header=FALSE))
  
  names(w1) <- names(w2) <- c("year",0:8)
  
  d <- reshape2::melt(w1,c("year"))
  names(d) <- c("year","age","w1")
  w2 <- reshape2::melt(w2,c("year"))
  d$w2 <- w2$value
  
  d$age <- as.integer(as.character(d$age))
  
  # Maturity
  p1 <- grep('#\tEnd\tof\tfile\tflag',txt) - 2
  p2 <- p1 + 1
  x <- str_clean(txt,p1,p2,header=T)
  x <- reshape2::melt(x)
  names(x) <- c("age","mat")
  x$age <- as.integer(as.character(x$age))
  d <- plyr::join(d,x,by="age")
  
  return(d)
}

#' @title read_sim_hake
#' 
#' @description Reads the simulation files for Hake.
#' 
#' @return Returns a list containing two data.frames, the raw data (raw) and
#' summarized data (qua). The latter contains the 0.05, 0.25, 0.50, 0.75 and
#' 0.95 percentiles as well as the mean.
#' 
#' @export
#' 
#' @param path The directory path where the Files are located files are located
#' 
read_sim_hake <- function(path) {
  
  # fingerplots
  filename <- dir(path,pattern="fingerplots",full.names = TRUE)
  p1 <- stringr::str_locate(filename,"fingerplots_")[,2] + 1
  p2 <- stringr::str_locate(filename,".out")[,1] - 1
  runname <- as.numeric(stringr::str_sub(filename,p1,p2))
  for (i in 1:length(filename)) {
    x <- read.table(filename[i])
    names(x) <- c("year","TAC","B_B1990","exp")
    n <- length(unique(x$year))
    iter <- nrow(x)/n
    x$iter <- rep(1:iter,each=n)
    x$run <- runname[i]
    if(i == 1) {
      raw <- x
    } else {
      raw <- rbind(raw,x)}
  }
  
  # Economic
  filename <- dir(path,pattern="Economic",full.names = TRUE)
  p1 <- stringr::str_locate(filename,"Economic_")[,2] + 1
  p2 <- stringr::str_locate(filename,".out")[,1] - 1
  runname <- as.numeric(stringr::str_sub(filename,p1,p2))
  for (i in 1:length(filename)) {
    x <- read.table(filename[i])
    names(x) <- c("year","TAC","value","employm","profit","vessel1","vessel2")
    n <- length(unique(x$year))
    iter <- nrow(x)/n
    x$iter <- rep(1:iter,each=n)
    x$run <- runname[i]
    if(i == 1) {
      raw2 <- x
    } else {
      raw2 <- rbind(raw2,x)}
  }
  raw <- plyr::join(raw,raw2[,-2],by=c("year","iter","run"))
  
  ## Depletion
  #filename <- dir(path,pattern="Depletion",full.names = TRUE)
  #p1 <- stringr::str_locate(filename,"Depletion_")[,2] + 1
  #p2 <- stringr::str_locate(filename,".out")[,1] - 1
  #runname <- as.numeric(stringr::str_sub(filename,p1,p2))
  #for (i in 1:length(filename)) {
  #  x <- read.table(filename[i])
  #  names(x) <- c("year","TAC","ssb","ssb_virgin")
  #  n <- length(unique(x$year))
  #  iter <- nrow(x)/n
  #  x$iter <- rep(1:iter,each=n)
  #  x$run <- runname[i]
  #  if(i == 1) {
  #    raw2 <- x
  #  } else {
  #    raw2 <- rbind(raw2,x)}
  #}
  #raw <- plyr::join(raw,raw2[,-2],by=c("year","iter","run"))
    
  q <- reshape2::melt(raw,c("year","iter","run"))
  q <- plyr::ddply(q,c("year","run","variable"),summarize,
             q05=quantile(value,0.05,na.rm=T),
             q25=quantile(value,0.25,na.rm=T),
             q50=quantile(value,0.50,na.rm=T),
             q75=quantile(value,0.75,na.rm=T),
             q95=quantile(value,0.95,na.rm=T),
             m=mean(value,na.rm=T))
  return(list(raw=raw,qua=q))
}

################################################################################
## Incomplete stuff below
#read_std_hake <- function(filename,runname) 
#{
#  d <- read.table(filename,header = TRUE,sep="")
#  return(d)
#}

##tBy/tB1 tBlastyear/B1990 AIKE spwnlastyear/spwnlast presentvalue
#read_stats_hake <- function(filename) 
#{
  #filename <- "ass/Results/stats.out"
#  txt <- readLines(filename)
#  txt <- str_trim(txt)
#  return(txt)
#}

#read_depletion_hake <- function(filename) {
#  d <- read.table(filename)
#  names(d) <- c("year","TAC","ssb","ssb_virgin")
#  n <- length(unique(d$year))
#  iter <- nrow(d)/n
#  d$iter <- rep(1:iter,each=n)
#  return(d)
#}

#read_depletion1990_hake <- function(filename) {
#  d <- read.table(filename)
#  names(d) <- c("year","TAC","bio_bio1990")
#  n <- length(unique(d$year))
#  iter <- nrow(d)/n
#  d$iter <- rep(1:iter,each=n)
#  return(d)
#}