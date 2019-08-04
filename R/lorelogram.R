#' Calculate pairwise log-odds ratios at intervals of increasing length
#'
#' This function estimates pairwise log-odds ratios in binary data at intervals of increasing length between subsequent sampling replicates. Plotting of the estimates versus lags (see the \code{\link{lor_plot}} function) provides a graphical description of how correlation between outcomes x-lag apart changes at the increase of the distance (in space or time) between the sampling replicates.
#'
#' @param data data.frame containing binary data in the wide format, with rows representing sampling units (e.g., camera trap sites or transects) and columns representing repeated samplings (e.g., temporal occasions or spatial replicates). See Details section below.
#' @param data_format character. Data can be provided in "wide" (default) or "long" format. See details for more information.
#' @param max_lag numeric. The maximum spatial or temporal lag between two sampling occasions at the same sampling units (default: 30) that should be considered.
#' @param lor_type character. Lorelogram can either be estimate using an "empirical" (default) or "model-based" approach.
#' @param id_rand_eff logical. Conditional on lor_type = "model-based". Does the model has to include sampling unit ID as random effect?
#' @param bin_width numeric. Number of lag that should be included in each bin. The default (=1) represents no binning.
#' @param plot_LOR logical. Create a plot of the results (default: TRUE)?
#' @param write_csv logical. Should the output be saved as a .csv (default: FALSE)?
#' @param outDir character. Directory into which .csv and plot file are saved.
#' @return A data.frame containing estimates of pairwise log-odds ratios and associated 95\% confidence intervals for each lag between 1 and \code{max_lag} is returned.
#' @details \code{\link{lorelogram}} can handle NAs in \code{data}. \code{data} should resemble a binary detection/nondetection history matrix such that provided as an output by the function \code{\link[camtrapR:detectionHistory]{camtrapR::detectionHistory}}. The first column should contain unitID (a unique identifier for each sampling units); each column from the second to the last should contain the binary data and should follow the spatial or temporal order in which the data were collected (e.g., second, third, and fourth columns should contain data from the first, second, third sampling replicates and so on). If unequal intervals are present in the data, fill gaps with columns of all NAs. Spatio- or temporal- difference between two subsequent columns (e.g. second and third) corresponds to 1-unit lag.
#'
#' @examples
#' data(GrayFox_Hour)
#' lorelogram(GrayFox_Hour, max_lag = 60)
#'
#'
#'
#' @importFrom magrittr %>%
#' @export
lorelogram <- function(data, data_format = "wide", max_lag = 30, lor_type = "empirical", id_rand_eff = FALSE, bin_width = 1, plot_LOR = TRUE, write_csv = FALSE, outDir = "") {

  wd0 <- getwd()
  on.exit(setwd(wd0))

  if (class(data)!="data.frame") {
    stop("Data should be a data.frame",
         call. = FALSE)
  }

  if (bin_width < 1) {
    stop("bin_width should be equal to 1 or higher",
         call. = FALSE)
  }


  if (data_format == "wide") {
  # Remove rows (=cameras) with no detection and prepare data
  y <- data[rowSums(data[,2:ncol(data)], na.rm = TRUE) > 0,]
  y <- droplevels(y)
  colnames(y) <- c("id", paste0("R", seq(1, J, 1), sep="")) # rename columns

  #y <- dplyr::rename(y, id=names(y[1]))

  # #### Organize data: from wide format to long format
  #+ organize1
  # Organize data
  y2 <- tidyr::gather(y, time, value, -id)
  y3 <- dplyr::mutate(y2, time=as.numeric(substr(time,2,10)))
  names(y3)[c(1,3)]<-c("id","y")
  y3$y <- as.numeric(as.character(y3$y))
  #head(y3)
  rm(list=c("y2"))

  # Define max number of reps
  max_reps <- ncol(y)
  } # close wide format

  if (data_format == "long") {
    y <- dplyr::rename(data, id = names(data[1]), time = names(data[2]), y = names(data[3]))
    # Remove rows (=cameras) with no detection and prepare data
    y3 <- as.data.frame(y %>% dplyr::group_by(id) %>% dplyr::filter(sum(y, na.rm = TRUE) > 0) %>% droplevels())
    y3[2:3] <- lapply(y3[2:3], as.numeric)

    # Define max number of reps
    max_reps <- max(y3$time, na.rm = TRUE)

  } # close long format


  # Determine all combinations of (current time, future times) for a sampling site up to max_lag.
  #
  # - V1 of x_cmb determines time point 1
  # - V2 of x_cmb determines time point 2
  #+ all_combinations

  if (max_reps>24*60*14) { #2-week data at minute time-interval
    # Create all combinations of minute-occasion up to 2 weeks of data
    x_cmb <- as.data.frame(arrangements::combinations(n=24*60*14, k=2, replace = FALSE))
    #head(x_cmb)

    # Now, calculate the time differences for all of these combinations and get rid of rows that include combinations of times where the difference in time <= max_lag.
    x_cmb <- x_cmb %>%
      dplyr::mutate(diff_time= x_cmb[,2] - x_cmb[,1]) %>% # column with time difference
      dplyr::filter(diff_time >= 0 & diff_time <= max_lag) %>% # filter combinations negative or too far apart
      dplyr::select(V1, V2)

    # Replicate to include more than 2 weeks of by-minute data
    x_cmb2 <- x_cmb
    for (i in 1:ceiling(max_reps/(24*60*14))){
      x_cmb_temp <- x_cmb + i*(24*60*14-max_lag)
      x_cmb2 <- rbind(x_cmb2,x_cmb_temp) %>% dplyr::filter(V1 <= max_reps)
    }
    x_cmb <- x_cmb2
    rm(list=c("x_cmb2", "x_cmb_temp", "data"))

    # Remove duplicates and interval larger than time of the last occasion
    x_cmb <- x_cmb %>% dplyr::distinct() %>% dplyr::filter(V1 <= max_reps & V2 < max_reps) %>% dplyr::mutate(diff_time = V2 - V1)
  } else {
    x_cmb <- as.data.frame(arrangements::combinations(n=max_reps, k=2, replace = FALSE))
    #head(x_cmb)

    # Now, calculate the time differences for all of these combinations and get rid of rows that include combinations of times where the difference in time <= max_lag.
    x_cmb <- x_cmb %>%
      dplyr::mutate(diff_time= x_cmb[,2] - x_cmb[,1]) %>% # column with time difference
      dplyr::filter(diff_time >= 0 & diff_time <= max_lag) %>% # filter combinations negative or too far apart
      dplyr::select(V1, V2)

    # Remove duplicates and interval larger than time last occasion
    x_cmb <- x_cmb %>% dplyr::distinct() %>% dplyr::filter(V1 <= max_reps & V2 < max_reps) %>% dplyr::mutate(diff_time = V2 - V1)
  }


  # #### Organize data: create pairwise detection histories
  # Create nested data frame with each row representing the data from a different time point.  This will make it easy to select data from all clusters where observations differ by a specific time lag.
  #+ organize2
  y.nest<-y3 %>% tidyr::nest(-time)
  #y.nest # nested data frame
  #as.data.frame(y.nest[1,]$data )[1:5,] # first 5 observations from first time point

  n <- ceiling(nrow(x_cmb)/1000000)
  a <- split(x_cmb, sort((1:nrow(x_cmb)) %% n))

  myfunct <- function(x) {
    # Sample first times and second times according to the first and second column of x_cmb, and include all clusters that have data from those time periods.  Then unnest to get one data set in long format.
    y.t1<- y.nest[a[[x]][,1],] %>% tidyr::unnest()
    y.t2<- y.nest[a[[x]][,2],] %>% tidyr::unnest()

    # Rename variables and cbind together
    names(y.t1)<-c("time1", "id", "y1")
    names(y.t2)<-c("time2", "id", "y2")
    Z<-cbind(y.t1[,c(2,1,3)], y.t2[,c(1,3)])
    Z<-Z %>% dplyr::mutate(time_diff=time2-time1) %>%
      dplyr::filter(!is.na(y1) & !is.na(y2))
    rm(list=c("y.t1","y.t2"))
    Z
  }

  b <- lapply(1:n, myfunct)
  Z <- as.data.frame(data.table::rbindlist(b, fill = TRUE)) # from list to data.frame

# #### Organize data: Compile counts of pairwise 11, 01, 01, 00 for each time interval
#+ organize3

if (lor_type == "empirical"){
  freq <- Z %>%
    dplyr::group_by(time_diff, y1, y2) %>%
    dplyr::summarise(count=n()) %>%
    tidyr::unite(col = y1y2, y1, y2, sep="", remove=FALSE) %>%
    dplyr::select(-y2) %>%  #remove info of second obs in time
    tidyr::spread(key=y1y2, value=count, fill=0)

    if (bin_width == 1){
      freq_tot <- freq %>%
        dplyr::group_by(time_diff, y1) %>%
        dplyr::summarise_all(dplyr::funs(sum))
      freq_tot
    }

    if (bin_width > 1) {
      # function from https://www.r-bloggers.com/finding-the-midpoint-when-creating-intervals/
      midpoints <- function(x, dp=2){
        lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
        upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
        return(round(lower+(upper-lower)/2, dp))
      }

      # build intervals
      brks = seq(min(freq$time_diff, na.rm = TRUE), max(freq$time_diff, na.rm = TRUE)+bin_width, bin_width)

      freq_tot <- freq %>%
        dplyr::mutate(bin = cut(time_diff, brks, include.lowest = TRUE),
                   Lag_midpoint = midpoints(bin)) %>%
        dplyr::select(-bin) %>%
        dplyr::group_by(Lag_midpoint) %>%
        dplyr::summarise_all(dplyr::funs(sum)) %>%
        dplyr::select(-time_diff) %>%
        dplyr::rename(., time_diff = Lag_midpoint)
      freq_tot
    }


# #### Calculate log odds ratios (Empirical lorelogram)
#+ Log_calculation
ORs <- freq_tot %>% dplyr::group_by(time_diff) %>%
  dplyr::summarize(or = sum(`11`)*sum(`00`)/(sum(`10`)*sum(`01`)),
                   se = sqrt(1/sum(`11`)+1/sum(`00`)+1/sum(`10`)+1/sum(`01`)))
ORs_all_minute <- ORs %>% dplyr::mutate(LORs=log(or),
                                        U_95_CI = log(or)+1.96*se,
                                        L_95_CI = log(or)-1.96*se)
ORs_all_minute[!is.na(ORs_all_minute$or) & ORs_all_minute$or==0,4:6] <- NA #replace values when lor values were undefined
LORs <- ORs_all_minute %>% dplyr::rename(., Lag=time_diff) %>% dplyr::select(c(Lag, LORs, U_95_CI, L_95_CI))
LORs
} # close lor_type empirical

if (lor_type == "model-based"){
  # organize the data
  if (id_rand_eff == FALSE){
    Z_par <- Z %>%
      dplyr::select(id, y1, y2, time_diff) %>%
      tidyr::unite(col = y2y1, y2, y1, sep="", remove=FALSE) %>%
      dplyr::group_by(time_diff, y1, y2y1) %>% #GROUP ID REMOVED
      dplyr::summarise(count = n()) %>%
      tidyr::spread(key = y2y1, value = count, fill=NA, sep = "_n")

    if (bin_width == 1){
      freq_tot <- Z_par %>%
        dplyr::mutate(success = ifelse(y1==1, y2y1_n11, y2y1_n10),                                      failure = ifelse(y1==1, y2y1_n01, y2y1_n00))
      freq_tot
    }

    if (bin_width > 1) {
      # function from https://www.r-bloggers.com/finding-the-midpoint-when-creating-intervals/
      midpoints <- function(x, dp=2){
        lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
        upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
        return(round(lower+(upper-lower)/2, dp))
      }
      # build intervals
      brks = seq(min(Z_par$time_diff, na.rm = TRUE), max(Z_par$time_diff, na.rm = TRUE)+bin_width, bin_width)
      freq_tot <- Z_par %>%
        dplyr::mutate(bin = cut(time_diff, brks, include.lowest = TRUE),
                      Lag_midpoint = midpoints(bin)) %>%
        dplyr::select(-bin) %>%
        dplyr::group_by(Lag_midpoint, y1) %>%
        dplyr::summarise_all(dplyr::funs(sum)) %>%
        dplyr::select(-time_diff) %>%
        dplyr::rename(., time_diff = Lag_midpoint) %>%
        dplyr::mutate(success = ifelse(y1==1, y2y1_n11, y2y1_n10),
                      failure = ifelse(y1==1, y2y1_n01, y2y1_n00))
      freq_tot
      }
    # run model and extract coefficients
    mod <- glmmTMB::glmmTMB(cbind(success, failure) ~ -1 + y1:as.factor(time_diff) + as.factor(time_diff), data = freq_tot, family = "binomial") # excluding random intercept
    LORs <- as.data.frame(cbind(Lag=unique(freq_tot$time_diff), confint(mod)[(length(unique(freq_tot$time_diff))+1):(nrow(confint(mod))),]))
    LORs <- dplyr::rename(LORs, L_95_CI=names(LORs[2]), U_95_CI=names(LORs[3]), LORs=names(LORs[4]))
    LORs
  } #close random effect false

  if (id_rand_eff == TRUE){
    Z_par <- Z %>%
      dplyr::select(id, y1, y2, time_diff) %>%
      tidyr::unite(col = y2y1, y2, y1, sep="", remove=FALSE) %>%
      dplyr::group_by(id, time_diff, y1, y2y1) %>% #GROUPED BY SAMPLING UNITS ID
      dplyr::summarise(count = n()) %>%
      tidyr::spread(key = y2y1, value = count, fill=NA, sep = "_n")

    if (bin_width == 1){
      freq_tot <- Z_par %>%
        dplyr::mutate(success = ifelse(y1==1, y2y1_n11, y2y1_n10),
                      failure = ifelse(y1==1, y2y1_n01, y2y1_n00))
      freq_tot
    }

    if (bin_width > 1) {
      # function from https://www.r-bloggers.com/finding-the-midpoint-when-creating-intervals/
      midpoints <- function(x, dp=2){
        lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
        upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
        return(round(lower+(upper-lower)/2, dp))
      }
      # build intervals
      brks = seq(min(Z_par$time_diff, na.rm = TRUE), max(Z_par$time_diff, na.rm = TRUE)+bin_width, bin_width)
      freq_tot <- Z_par %>%
        dplyr::mutate(bin = cut(time_diff, brks, include.lowest = TRUE),
                      Lag_midpoint = midpoints(bin)) %>%
        dplyr::select(-bin) %>%
        dplyr::group_by(id, Lag_midpoint, y1) %>%
        dplyr::summarise_all(dplyr::funs(sum)) %>%
        dplyr::select(-time_diff) %>%
        dplyr::rename(., time_diff = Lag_midpoint) %>%
        dplyr::mutate(success = ifelse(y1==1, y2y1_n11, y2y1_n10),
                      failure = ifelse(y1==1, y2y1_n01, y2y1_n00))
      freq_tot
    }

    # run model and extract coefficients
    mod <- glmmTMB::glmmTMB(cbind(success, failure) ~ -1 + y1:as.factor(time_diff) + as.factor(time_diff) + (1|id), data = freq_tot, family = "binomial") # including random intercept
    LORs <- as.data.frame(cbind(Lag=unique(freq_tot$time_diff), confint(mod)[(length(unique(freq_tot$time_diff))+1):(nrow(confint(mod))-1),]))
    LORs <- dplyr::rename(LORs, L_95_CI=names(LORs[2]), U_95_CI=names(LORs[3]), LORs=names(LORs[4]))
    LORs
  } # close random effect true

  LORs
} # close lor_type model-based

if (write_csv == TRUE) { # save results

  if(outDir != "" & !file.exists(outDir)){stop("outDir does not exist")}

  if (outDir == "") {
    filename <- paste("Log_odds_ratio_MaxLag_", as.character(max_lag),".csv", sep="")
    } else {
    filename <- paste(as.character(outDir), "/Log_odds_ratio_MaxLag_", as.character(max_lag),".csv", sep="")
    }
  write.csv(x = LORs, file = filename)
}

if (plot_LOR == TRUE) { #plot lor values
  plot_LOR_all <- lor_plot(LORs)
  print(plot_LOR_all)
}

LORs

}
