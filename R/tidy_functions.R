#' @title str_clean
#' 
#' @description A helper function
#' 
#' @param txt A text vector
#' @param p1 Numerical value specifying 1st position
#' @param p2 Numerical value specifying last position
#' @param Replace A text string
#' @param header Boolean
#' @param sep A text string
#' @param ... Additional arguments


str_clean <- function(txt,p1,p2,Replace='',header=FALSE,sep='', ...) {
  txt <- stringr::str_trim(txt)
  txt <- txt[(p1):(p2)]
  txt <- stringr::str_replace(txt,'#',Replace)
  txt <- str_trim_tab(txt)
  tmpfile <- tempfile()
  writeLines(txt,tmpfile)
  df <- read.table(tmpfile,header=header,sep=sep,fill=TRUE)
  names(df) <- stringr::str_replace(names(df),'X','')
  return(df)
}

#' @title str_trim_tab
#' 
#' @title Trim tabs from start and end of string
#' 
#' @param txt input character vector
#' @param side side on which tab is removed (left,right,both)
#' 
str_trim_tab <- function (txt, side = "both") {
  #string <- check_string(string)
  stopifnot(length(side) == 1)
  side <- match.arg(side, c("left", "right", "both"))
  pattern <- switch(side,
                    left = "^\\t+",
                    right = "\\t+$", 
                    both = "^\\t+|\\t+$")
  stringr::str_replace_all(txt, pattern, "")
}