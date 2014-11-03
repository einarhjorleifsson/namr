#' @title yield per recruit by age
#' 
#' @description Returns a detailed (by age) yield and ssb composition. 
#' Summarization over age is provided by the function \code{ypr}
#' 
#' @export
#' 
#' @param x \code{data.frame}
#' @param fmort A numerical value specifying fishing mortality
#' @param a1 Value containing lower reference age
#' @param a2 Value containing upper reference age
#' @param recruitment An number specifying number of recruits. If one (default)
#' returns the conventional yield and ssb for one recruit.
#' @param plusgroup A boolean. If TRUE (default) the last age group is treated
#' as a plus group.
#' @param pM Value containing the proportion of M prior to spawning
#' @param pF Value containing the proportion of F prior to spawning
#' 
#' @return returns \code{list} with catch and ssb for each age
#'
ypra <- function(x,
                 fmort=0,
                 a1,
                 a2,
                 plusgroup=TRUE,
                 recruitment=1,
                 pM=0,
                 pF=0)
{
  
  if(is.null(x$pM)) x$pM <- pM
  if(is.null(x$pF)) x$pF <- pF
  
  # renormalise selection in case sel is fishing mortality
  x$sel <- x$sel/mean(x$sel[x$age %in% a1:a2])
  
  x$f <- x$sel * fmort
  #x$z <- x$f + x$m
  x$c <- x$nSSB  <- x$n <- NA
  
  x$n[1] <- recruitment
  x$nSSB[1] <- x$n[1] * exp(-(x$m[1]*x$pF[1] + x$f[1]*x$pM[1]))
  x$c[1] <- x$n[1] * exp(-0.5 *x$m[1]) *  x$f[1]

  n <- nrow(x)
  
  for(i in 1:(n-1)) {
    x$n[i+1] <- (x$n[i] * exp(-0.5 * x$m[i]) - x$c[i]) * exp(-0.5 * x$m[i])
    x$c[i+1] <- x$n[i+1] * exp(-0.5 * x$m[i+1]) * x$f[i+1]
    x$nSSB[i+1] <- x$n[i+1] * exp(-(x$m[i+1]*x$pM[i+1] + x$f[i+1]*x$pF[i+1]))
  }
  
  if(plusgroup) {
    nplus <- (x$n[n] * exp(-x$m[n]/2) - x$c[n]) * exp(-x$m[n]/2)
    q <- nplus/x$n[n]
    x$n[n]  <- x$n[n]/(1-q)
    x$nSSB[n] <- x$n[n] * exp(-(x$m[i+1]*x$pM[i+1] + x$f[i+1]*x$pF[i+1]))
    x$c[n] <- x$c[n]/(1-q)
  }

  x$bio <- x$n * x$sW
  x$ssb <- x$nSSB * x$sW * x$mat
  x$exp <- x$n * exp(-0.5 * x$m) * x$cW * x$sel
  x$yield <- x$c * x$cW
  return(x)
}

#' @title yield per recruit
#' 
#' @description Calculates yield and ssb per recruit. 
#' Summarization over age is provided by the function \code{ypra}
#' 
#' @export
#' 
#' @param df \code{data.frame}
#' @param fscan A numerical vector specifying fishing mortality to scan over
#' @param a1 Value containing lower reference age
#' @param a2 Value containing upper reference age
#' @param plusgroup A boolean. If TRUE (default) the last age group is treated
#' as a plus group.
#' @param recruitment An number specifying number of recruits. If one (default)
#' returns the conventional yield and ssb for one recruit.
#' @param pM Value containing the proportion of M prior to spawning
#' @param pF Value containing the proportion of F prior to spawning
#' 
#' @return returns \code{list} with catch and ssb for each age
#'
ypr <- function(df,
                fscan=0,
                a1,
                a2,
                plusgroup=TRUE,
                recruitment=1,
                pM=0,
                pF=0)
{
  
  # dummy stuff
  #yield <- ssb <- NULL
  
  if(is.null(df$pM)) df$pM <- pM
  if(is.null(df$pF)) df$pF <- pF
  
  # renormalise selection in case sel is fishing mortality
  df$sel <- df$sel/mean(df$sel[df$age %in% a1:a2])
  
  n <- length(fscan)  
  
  yprbyage <- data.frame(f=rep(fscan,each=nrow(df)),
                         age=rep(df$age,n),
                         m=rep(df$m,n),
                         sel=rep(df$sel,n),
                         pF=rep(df$pF,n),
                         pM=rep(df$pM,n),
                         mat=rep(df$mat,n),
                         cW=rep(df$cW,n),
                         sW=rep(df$sW,n))
  yprbyage$yield <- yprbyage$ssb <- NA
  
  for (i in 1:length(fscan)) {
    j <- yprbyage$f == fscan[i]
    tmp <- ypra(yprbyage[j,],fmort=fscan[i],a1=a1,a2=a2,plusgroup=plusgroup,recruitment=recruitment)
    yprbyage$yield[j] <- tmp$yield
    yprbyage$ssb[j]   <- tmp$ssb
  }
  
  res <- plyr::ddply(yprbyage,c("f"),summarise,yield=sum(yield),ssb=sum(ssb))
  
  x <- spline(res$f,res$yield,n=10*nrow(res))
  maxy <- max(x$y)
  fmax <- x$x[x$y==max(x$y)]
  d <- diff(x$y)
  d <- d/d[1]
  d1 <- d[1:(length(d)-1)]
  d2 <- d[2:length(d)]
  i <- 1:length(d1)
  i <- i[d1 > 0.1 & d2 < 0.1]
  f01 <- x$x[i]
  
  x <- spline(res$f,res$ssb,n=10*nrow(res))
  x1  <- abs(x$y-0.35*max(res$ssb,na.rm=T)) # check why returns one NaN
  ssb35 <- x$x[x1==min(x1)]
  
  refpts <- rep(0,4)
  names(refpts) <- c("f01","fmax","maxy","ssb35")
  refpts["f01"] <- f01
  refpts["fmax"] <- fmax
  refpts["maxy"] <- maxy
  refpts["ssb35"] <- ssb35
  refpts <- reshape2::melt(refpts)
  refpts$variable <- rownames(refpts)
  rownames(refpts) <- NULL
  refpts <- refpts[,c("variable","value")]
  
  return(list(ypr=res,yprba=yprbyage,refpts=refpts))
}



