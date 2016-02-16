# To Do
# * Calculate oC and pC as a sum for the whole fleet and put into rbya


#' @title Read Horsemackerel run
#' 
#' @description XXX
#' 
#' @param filename Normally ends on .rep
#' 
#' @export

read_hmac <- function(filename) {
 
  rbx <- read_rbx_hmac(filename)
  fn  <- stringr::str_replace(filename, ".rep", "data.dat")
  ibx <- read_ibx_hmac(fn)
  
  # By year and age
  rbya <- 
    rbx$rbya %>% 
    dplyr::left_join(ibx$ibya, by = c("year", "age")) %>% 
    dplyr::select(year, age, cW1, cW2, n, f) %>% 
    dplyr::left_join(ibx$iba, by = "age") # include the maturity
  rbya$m <- 0.45   # Natural mortality
  # By year and fleet
  rbyf <- 
    rbx$rbyf %>%
    dplyr::left_join(ibx$ibyf, by = c("year", "fleet"))

  # By year
  #  sum landings by fleet
  rby <- 
    rbyf %>% 
    dplyr::select(year, landings) %>% 
    dplyr::group_by(year) %>%
    dplyr::summarise(landings = sum(landings, na.rm = TRUE)) %>% 
    dplyr::ungroup() %>% 
    dplyr::right_join(rbx$rby, by = "year")
    
  return(list(rby = rby, rbya = rbya, rbyf = rbyf, rbaf = rbx$rbaf, rbyaf = rbx$rbyaf))
  
}


#' Read Horse Mackerel configure file
#' 
#' @description XXXX
#' 
#' @param filename XXX
#' 
#' @export
#' 
read_conf_hmac <- function(filename) {
  
  cntr <- read.table(filename,
                     col.names = c("value","cntr"),
                     skip = 3, 
                     comment.char = "") %>% 
    dplyr::mutate(cntr = stringr::str_replace(cntr, "#", ""))
  ctr <- as.list(cntr$value)
  names(ctr) <- cntr$cntr
  
  return(ctr)
}

#' Title
#'
#' @param filename XXX
#'
#' @export

read_ibx_hmac <- function(filename) {

  # ----------------------------------------------------------------------------
  # Constants
  AGES = c(0:7)
  FLEETS=c('Midwater','Pelagic','Bulgaria','Poland','Rumania','USSR','Survey')
  
  # ----------------------------------------------------------------------------
  # Controls
  ctr.file <- stringr::str_replace(filename, "data", "")
  ctr  <- read_conf_hmac(ctr.file)
  YRS1 <- c(ctr$first_yr:ctr$last_yr)
  YRS2 <- c(YRS1, c(ctr$last_yr + c(1:(ctr$nfut_yrs+1))))
  
  txt <- readLines(filename)
  txt <- str_trim_tab(txt) # not only trims tabs at both ends
  
  # ----------------------------------------------------------------------------
  # By age

  # Maturity at age:
  p1 <- grep('#maturity',txt)
  p2 <- p1 + 2
  d <- txt[p2]
  iba <-  
    as.data.frame(d) %>% 
    tidyr::separate(d, AGES, sep = "\t", convert = TRUE) %>% 
    tidyr::gather(age, mat) %>% 
    mutate(age = as.integer(age))
  
  # ----------------------------------------------------------------------------
  # By year and fleet

  # Landings
  p1 <- grep('#\tMidwater',txt)
  p2 <- grep('#Note:',txt)
  d <- txt[(p1+2):(p2-2)]
  ibyf <-
    as.data.frame(d) %>% 
    tidyr::separate(d, YRS1, sep = "\t", convert = TRUE) %>% 
    dplyr::mutate(fleet = FLEETS[-7]) %>% 
    tidyr::gather(year, landings, -fleet, convert = TRUE)
  
  # ----------------------------------------------------------------------------
  # By year and age and fleet

  # Proportion of numbers caught by age by midwater fleet (Namibia)
  p1 <- grep('#\tProportion\tof\tNumbers\tcaught\tby\tage\tby\tMidwater',txt)
  p2 <- grep('#Proportion\tof\tNumbers\tcaught\tby\tthe\tpelagic',txt)
  d <- txt[(p1+2):(p2-2)]
  ibyaf <-
    as.data.frame(d) %>% 
    tidyr::separate(d, c("year",AGES), sep = "\t", convert=TRUE) %>% 
    tidyr::gather(age, oC, -year, convert = TRUE) %>% 
    dplyr::mutate(fleet = "Midwater")

  # Proportion of numbers caught by age by pelagic fleet (Namibia)
  p1 <- grep('#Proportion\tof\tNumbers\tcaught\tby\tthe\tpelagic',txt)
  p2 <- grep('#Bulgaria\tcatch-at-age',txt)
  d <- txt[(p1+2):(p2-2)]
  ibyaf <-
    as.data.frame(d) %>% 
    tidyr::separate(d, c("year",AGES), sep = "\t", convert=TRUE) %>% 
    tidyr::gather(age, oC, -year, convert = TRUE) %>% 
    dplyr::mutate(fleet = "Pelagic") %>% 
    dplyr::bind_rows(ibyaf)
  
  # Proportion of numbers caught by age by Bulgaria
  p1 <- grep('#Bulgaria\tcatch-at-age',txt)
  p2 <- grep('#Poland\tcatch-at-age',txt)
  d <- txt[(p1+2):(p2-2)]
  ibyaf <-
    as.data.frame(d) %>% 
    tidyr::separate(d, c("year",AGES), sep = "\t", convert=TRUE) %>% 
    tidyr::gather(age, oC, -year, convert = TRUE) %>% 
    dplyr::mutate(fleet = "Bulgaria") %>% 
    dplyr::bind_rows(ibyaf)  

  # Proportion of numbers caught by age by Poland
  p1 <- grep('#Poland\tcatch-at-age',txt)
  p2 <- grep('#Romania\tcatch-at-age',txt)
  d <- txt[(p1+2):(p2-2)]
  ibyaf <-
    as.data.frame(d) %>% 
    tidyr::separate(d, c("year",AGES), sep = "\t", convert=TRUE) %>% 
    tidyr::gather(age, oC, -year, convert = TRUE) %>% 
    dplyr::mutate(fleet = "Poland") %>% 
    dplyr::bind_rows(ibyaf)
  
  # Proportion of numbers caught by age by Romania
  p1 <- grep('#Romania\tcatch-at-age',txt)
  p2 <- grep('#USSR\tcatch-at-age',txt)
  d <- txt[(p1+2):(p2-2)]
  ibyaf <-
    as.data.frame(d) %>% 
    tidyr::separate(d, c("year",AGES), sep = "\t", convert=TRUE) %>% 
    tidyr::gather(age, oC, -year, convert = TRUE) %>% 
    dplyr::mutate(fleet = "Romania") %>% 
    dplyr::bind_rows(ibyaf)
  
  # Proportion of numbers caught by age by USSR
  p1 <- grep('#USSR\tcatch-at-age',txt)
  p2 <- grep('#\tProportion\tof\tNumbers\tcaught\tby\tsurveys\tby\tyear',txt)
  d <- txt[(p1+2):(p2-2)]
  ibyaf <-
    as.data.frame(d) %>% 
    tidyr::separate(d, c("year",AGES), sep = "\t", convert=TRUE) %>% 
    tidyr::gather(age, oC, -year, convert = TRUE) %>% 
    dplyr::mutate(fleet = "USSR") %>% 
    dplyr::bind_rows(ibyaf)
  
  # Proportion of numbers caught by age in surveys
  p1 <- grep('#\tProportion\tof\tNumbers\tcaught\tby\tsurveys\tby\tyear',txt)
  p2 <- grep('#\tCommercial\tCPUE:',txt)
  d <- txt[(p1+2):(p2-2)]
  ibyaf <-
    as.data.frame(d) %>% 
    tidyr::separate(d, c("year",AGES), sep = "\t", convert=TRUE) %>% 
    tidyr::gather(age, oC, -year, convert = TRUE) %>% 
    dplyr::mutate(fleet = "Survey") %>% 
    dplyr::bind_rows(ibyaf)
  
  # ----------------------------------------------------------------------------
  # By year and fleet 2

  # Commerical catch per unit effort by fleet
  p1 <- grep('#\tCommercial\tCPUE:',txt)
  p2 <- grep('#\tsurvey\tseries',txt)
  d <- txt[(p1+3):(p2-2)]
  d <- d[c(1,3,5,7,9)]
  rbyf <-
    as.data.frame(d) %>% 
    tidyr::separate(d, YRS1, sep = "\t", convert = TRUE) %>% 
    dplyr::mutate(fleet = c('Midwater','Bulgaria','Poland','Romania','USSR')) %>% 
    tidyr::gather(year, cpue, -fleet, convert = TRUE) %>% 
    dplyr::mutate(cpue = ifelse(cpue == 0, NA, cpue))
  
  # ----------------------------------------------------------------------------
  # Survey total biomass indices
  p1 <- grep('#\tsurvey\tseries',txt)
  p2 <- grep('#startWeight',txt)
  d <- txt[(p1+2):(p2-2)]
  x <-
    as.data.frame(d) %>% 
    tidyr::separate(d, YRS1, sep = "\t", convert = TRUE)
  ibySurvey <-
    data.frame(year = as.integer(names(x)),
               oU = as.numeric(x[1,]),
               cvU = as.numeric(x[2,]))
  
  # ----------------------------------------------------------------------------
  # by year and age
  
  # The weight in start of year
  p1 <- p2 <- grep('#startWeight',txt)
  p2 <- grep('#midWeight',txt)
  d <- txt[(p1+2):(p2-1)]
  ibya <- 
    as.data.frame(d) %>% 
    tidyr::separate(d, c("year", AGES), sep = "\t", convert = TRUE) %>% 
    tidyr::gather(age, cW1, -year, convert = TRUE)

  # The weight in middle of year
  p1 <- grep('#midWeight',txt)
  p2 <- grep('#========Selectivities',txt)
  d <- txt[(p1+2):(p2-1)]
  ibya <- 
    as.data.frame(d) %>% 
    tidyr::separate(d, c("year", AGES), sep = "\t", convert = TRUE) %>% 
    tidyr::gather(age, cW2, -year, convert = TRUE) %>% 
    dplyr::left_join(ibya, by = c("year", "age"))

  # ----------------------------------------------------------------------------
  # By age and fleet
  p1 <- p2
  p2 <- grep('#----------------------------',txt)[2]
  d <- txt[(p1+2):(p2-1)]
  ibaf <- 
    as.data.frame(d) %>% 
    tidyr::separate(d, AGES, sep = "\t", convert = TRUE) %>% 
    dplyr::mutate(fleet = FLEETS) %>% 
    tidyr::gather(sel, age, -fleet, convert = TRUE)

  results <- list(iba   = dplyr::tbl_df(iba),
                  ibaf  = dplyr::tbl_df(ibaf),
                  ibya  = dplyr::tbl_df(ibya),
                  ibyaf = dplyr::tbl_df(ibyaf),
                  ibyf  = dplyr::tbl_df(ibyf),
                  ibySurvey = dplyr::tbl_df(ibySurvey))
  return(results)
}


#' @title Reads the ADMB report (result) file from the Namibian horse mackerel assessment
#'
#' @description XXX
#'
#' @export
#'
#' @param filename A character string, including path (directory) to the report
#' file. The report file has normally the ending .rep.

read_rbx_hmac <- function(filename) {
  
  # ----------------------------------------------------------------------------
  # Constants
  AGES = c(0:7)
  FLEETS=c('Midwater','Pelagic','Bulgaria','Poland','Rumania','USSR','Survey')
  
  # ----------------------------------------------------------------------------
  # Controls
  ctr.file <- stringr::str_replace(filename, ".rep", ".dat")
  ctr  <- read_conf_hmac(ctr.file)
  YRS1 <- c(ctr$first_yr:ctr$last_yr)
  YRS2 <- c(YRS1, c(ctr$last_yr + c(1:(ctr$nfut_yrs+1))))
  
  # ----------------------------------------------------------------------------
  # Read file
  txt <- readLines(filename)
  txt <- str_trim_tab(txt)
  
  # ----------------------------------------------------------------------------
  # AGE STUFF
  
  # Virgin age structure
  p1 <- grep('Virgin Age Structure',txt) + 1
  d <- txt[(p1)]
  d <- str_trim_all(d)
  rba <- data.frame(age = AGES, N0 = as.numeric(stringr::str_split(d, " ")[[1]]))

  # Fleet selectivity
  p1 <- grep('Selection by fleet by age. Order of fleets',txt) + 2
  p2 <- p1 + 5
  d <- txt[p1:p2]
  # Survey selectivity
  p1 <- grep('Survey selection by age',txt) + 1
  x <- txt[p1:p1]
  # combine and generate dataframe
  d <- c(d,x)
  rbaf <- as.data.frame(d) %>% 
    dplyr::mutate(d = stringr::str_trim(d),
                  fleet = FLEETS) %>% 
    tidyr::separate(d, AGES, sep = " ", convert = TRUE) %>% 
    tidyr::gather(age, sel, -fleet, convert = TRUE)
  
  # ----------------------------------------------------------------------------
  # Results by year and age

  # Stock in numbers
  p1 <- grep('N-matrix',txt)
  p2 <- grep('Fishing mortality per year per fleet',txt)
  d <- txt[(p1+1):(p2-2)]
  rbya <- as.data.frame(d) %>% 
    dplyr::mutate(d = stringr::str_trim(d),
                  year = YRS2) %>% 
    tidyr::separate(d, AGES, sep = " ", convert = TRUE) %>% 
    tidyr::gather(age, n, -year, convert = TRUE)
  
  # ----------------------------------------------------------------------------
  # Results by year and fleet 1
  
  # Fishing mortality per year per fleet 
  p1 <- grep('Fishing mortality per year per fleet',txt)
  p2 <- grep('Year    TotalB  Spawning biomass  Depletion',txt)
  d <- txt[(p1+2):(p2-2)]
  rbyf <- as.data.frame(d) %>% 
    dplyr::mutate(d = stringr::str_trim(d),
                  fleet = FLEETS[1:6]) %>% 
    tidyr::separate(d, YRS2, sep = " ", convert = TRUE) %>% 
    tidyr::gather(year, f, -fleet, convert = TRUE) %>% 
    dplyr::tbl_df()
  
  # ----------------------------------------------------------------------------
  # Results by year 1
  
  # Biomass information
  p1 <- grep('Year    TotalB  Spawning biomass  Depletion',txt)
  p2 <- grep('RY  Yield',txt)
  d <- txt[(p1+1):(p2-2)]
  rby <- as.data.frame(d) %>% 
    dplyr::mutate(d = stringr::str_trim(d)) %>% 
    tidyr::separate(d, c("year", "tBio", "ssb", "depletionBio1961", "depletionBio1990"),
                    sep = " ",
                    convert = TRUE)

  # RY  Yield(1,Year)  Yield(2,Year)
  #   RY is most likely replacement yield
  #   Unsure what yield1 and yield2 means
  #  For the time being this is put into object rby
  p1 <- grep(' RY  Yield',txt)
  p2 <- grep('Year - Exploitable biomass',txt)
  d <- txt[(p1+1):(p2-1)]
  d <- as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "rY", "y1", "y2"), sep = " ", convert = TRUE)
  
  rby <-
    dplyr::left_join(rby, d, by = "year")
  
  # ----------------------------------------------------------------------------
  # Results by year and fleet 2

  # Exploitable biomass for each fleet
  p1 <- grep('Year - Exploitable biomass',txt)
  p2 <- grep('Year       FutCatch',txt)
  d <- txt[(p1+1):(p2-2)]
  rbyf <- as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", FLEETS[-7]), sep = " ", convert = TRUE) %>% 
    tidyr::gather(fleet, explBio, -year, convert = TRUE) %>% 
    dplyr::left_join(rbyf, by = c("year", "fleet"))
  
  # ----------------------------------------------------------------------------
  # Results by year 2

  # Observed and predicted survey estimates
  #   The survey cv is most likely an input value
  p1 <- grep('Year Survey Est Survey val',txt)
  p2 <- grep('Year Obs pred',txt)
  d <- txt[(p1+1):(p2-1)]
  rby <- as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "oU", "pU", "cvU"), sep = " ", convert = TRUE) %>% 
    dplyr::mutate(oU = ifelse(oU == 0, NA, oU),
                  cvU = ifelse(cvU == 0, NA, cvU)) %>% 
    dplyr::right_join(rby, by = "year")

  # ----------------------------------------------------------------------------
  # Results by year and fleet 2
  
  p1 <- grep('Year Obs pred',txt)
  p2 <- grep('Recruitment residuals',txt)[2]
  
  # CV of the cpue is the first line (data are not returned in the function output)
  d <- txt[p1+1]
  cvFleet <- as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, FLEETS[c(1,3:6)], sep = " ", convert = TRUE)
  
  # Catch per unit effort
  # Fleets: Midwater, Bulgaria, Poland, Rumania, USSR
  # For each fleet we have observed and then predicted (0 means NA)
  cn1 <- rep(c("o", "p"), 5)
  cn2 <- rep(FLEETS[c(1,3:6)],each=2)
  cn  <- paste0(cn1,cn2)
  d <- txt[(p1+2):(p2-1)]
  d <- as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d),
                  d = stringr::str_replace_all(d, " 0 "," NA ")) %>% 
    tidyr::separate(d, c("year", cn), sep = " ", convert = TRUE)
  # a) take the observed and gather
  obs <-
    d %>% 
    dplyr::select(year, starts_with("o")) %>% 
    tidyr::gather(fleet, oU, -year, convert = TRUE) %>% 
    dplyr::mutate(fleet = stringr::str_sub(fleet, 2))
  # a) take the predicted and gather
  pre <-
    d %>% 
    dplyr::select(year, starts_with("p")) %>% 
    tidyr::gather(fleet, pU, -year, convert = TRUE) %>% 
    dplyr::mutate(fleet = stringr::str_sub(fleet, 2))
  # c) join with existing rbyf
  rbyf <-
    rbyf %>% 
    dplyr::left_join(obs, by = c("year", "fleet")) %>% 
    dplyr::left_join(pre, by = c("year", "fleet"))
  
  # ----------------------------------------------------------------------------
  # Results by year 3
  # Recruitment residuals
  p1 <- grep('Recruitment residuals',txt)[2]
  p2 <- grep('Fit to CAA Midwater data',txt)
  d <- txt[(p1+1):(p2-1)]
  rby <- 
    as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "rRec"), sep = " ", convert = T) %>% 
    dplyr::right_join(rby, by = "year")
  
  # ----------------------------------------------------------------------------
  # Results by year and age and fleet - CATCH AT AGE
  
  # Midwater fleet
  p1 <- grep('Fit to CAA Midwater data',txt)
  p2 <- grep('Fit to CAA Pelagic data',txt)
  d <- txt[(p1+1):(p2-1)]
  rbyaf <-
    as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "age", "oC", "pC", "rC"), sep = " ", convert = TRUE) %>% 
    dplyr::mutate(fleet = "Midwater")
  
  # Catch at age - Pelagic data
  p1 <- grep('Fit to CAA Pelagic data',txt)
  p2 <- grep('Fit to CAA Bulgaria data',txt)
  d <- txt[(p1+1):(p2-1)]
  rbyaf <-
    as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "age", "oC", "pC", "rC"), sep = " ", convert = TRUE) %>% 
    dplyr::mutate(fleet = "Pelagic") %>% 
    dplyr::bind_rows(rbyaf)

  # Catch at age - Bulgaria data
  p1 <- grep('Fit to CAA Bulgaria data',txt)
  p2 <- grep('Fit to CAA Poland data',txt)
  d <- txt[(p1+1):(p2-1)]
  rbyaf <-
    as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "age", "oC", "pC", "rC"), sep = " ", convert = TRUE) %>% 
    dplyr::mutate(fleet = "Bulgaria") %>% 
    dplyr::bind_rows(rbyaf)
  
  # Catch at age - Poland data
  p1 <- grep('Fit to CAA Poland data',txt)
  p2 <- grep('Fit to CAA Romania data',txt)
  d <- txt[(p1+1):(p2-1)]
  rbyaf <-
    as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "age", "oC", "pC", "rC"), sep = " ", convert = TRUE) %>% 
    dplyr::mutate(fleet = "Poland") %>% 
    dplyr::bind_rows(rbyaf)
  
  # Catch at age -  Romania data
  p1 <- grep('Fit to CAA Romania data',txt)
  p2 <- grep('Fit to CAA USSR data',txt)
  d <- txt[(p1+1):(p2-1)]
  rbyaf <-
    as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "age", "oC", "pC", "rC"), sep = " ", convert = TRUE) %>% 
    dplyr::mutate(fleet = "Romania") %>% 
    dplyr::bind_rows(rbyaf)
  
  # Catch at age - USSR data
  p1 <- grep('Fit to CAA USSR data',txt)
  p2 <- grep('Fit to CAA survey data',txt)
  d <- txt[(p1+1):(p2-1)]
  rbyaf <-
    as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "age", "oC", "pC", "rC"), sep = " ", convert = TRUE) %>% 
    dplyr::mutate(fleet = "USSR") %>% 
    dplyr::bind_rows(rbyaf)
  
  # Catch at age - Survey data
  p1 <- grep('Fit to CAA survey data',txt)
  p2 <- grep('Selection by fleet by age',txt)
  d <- txt[(p1+1):(p2-2)]
  rbyaf <-
    as.data.frame(d) %>% 
    dplyr::mutate(d = str_trim_all(d)) %>% 
    tidyr::separate(d, c("year", "age", "oC", "pC", "rC"), sep = " ", convert = TRUE) %>% 
    dplyr::mutate(fleet = "Survey") %>% 
    dplyr::bind_rows(rbyaf)
  
  # ----------------------------------------------------------------------------
  # Calculate Fay for each fleet
  Fy <- 
    rbyf %>% 
    dplyr::select(fleet, year, f)
  Fa <-
    rbaf %>% 
    dplyr::filter(fleet %in% unique(Fy$fleet))
  Fay <- 
    dplyr::left_join(Fy, Fa, by = c("fleet")) %>% 
    dplyr::mutate(f = sel * f) %>% 
    dplyr::select(-sel)
  
  # EINAR: Something is rotten in the State of Denmark
  #   The observed catches are at maximum equal to 1
  #   The sum of those observed catches are always equal to 1
  #   But the fishing mortality goes all the way to age 7!!
  rbyaf <- 
    rbyaf %>% 
    dplyr::right_join(Fay, by = c("year", "age", "fleet"))
  
  # ----------------------------------------------------------------------------
  # Calculate the Fay, i.e. sum of fleet F for each year
  rbya <- 
    Fay %>% 
    dplyr::group_by(year, age) %>% 
    dplyr::summarise(f = sum(f)) %>% 
    dplyr::ungroup() %>% 
    dplyr::right_join(rbya, by = c("year", "age"))
  
  
  # ----------------------------------------------------------------------------
  # Put survey oU and pU with other fleets
  x <-
    rby %>% 
    dplyr::mutate(fleet = "Survey") %>%
    dplyr::select(year, fleet, oU, pU)
  rbyf <-
    dplyr::bind_rows(rbyf,x)
  
  rby <-
    rby %>% 
    dplyr::select(-oU, -pU)
  
  # ----------------------------------------------------------------------------
  # Add recruitment to rby
  rec <- 
    rbya %>% 
    dplyr::filter(age == 0) %>% 
    dplyr::select(year, rec = n)
  rby <-
    rby %>% 
    dplyr::left_join(rec, by = "year")
  # reorder
  rby <- 
    rby %>% 
    dplyr::select(year, rec, bio = tBio, ssb, rRec, dBio1961 = depletionBio1961,
                  dBio1990 = depletionBio1990, rY, y1, y2)
  
  # ----------------------------------------------------------------------------

  results <- list(rbaf   = dplyr::tbl_df(rbaf),
                  rbya  = dplyr::tbl_df(rbya),
                  rbyf  = dplyr::tbl_df(rbyf),
                  rby   = dplyr::tbl_df(rby),
                  rbyaf = dplyr::tbl_df(rbyaf))
  return(results)
}
