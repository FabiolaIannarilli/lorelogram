#' Calculate and plot pairwise log-odds ratios at intervals of increasing length
#'
#' This function estimates pairwise log-odds ratios in binary data at intervals of increasing length between subsequent sampling replicates. Plot of the estimates versus lags provides a graphical description of how correlation between outcomes x-lag apart changes at the increase of the distance (in space or time) between the sampling replicates.
#'
#' @param data data.frame containing binary data in the wide format, with rows representing sampling units (e.g., camera trap sites or transects) and columns representing repeated samplings (e.g., temporal occasions or spatial replicates).
#' @param max_lag numeric. The maximum spatial or temporal lag between two sampling occasions at the same sampling units (default: 30).
#' @param write_csv logical. Should output be saved as a .csv (default: TRUE)?
#' @param outDir character. Directory into which .csv and plot file are saved.
#' @param plot_LOR logical. Create a .jpg plot of the results (default: TRUE)?
#' @param plot_title character. Title of the plot (default: NULL)
#' @param plot_x_title character. Title of x-axis of the plot (default: Lag)
#' @return A data.frame containing estimates of pairwise log-odds ratios and associated 95\% confidence intervals for each lag between 1 and \code{max_lag} is returned.
#' @details \code{data} should resemble a binary detection/nondetection history matrix such that provided as an output by the function \code{\link[camtrapR:detectionHistory]{camtrapR::detectionHistory}}. \code{\link{lorelogram}} can handle NAs in \code{data}.
#'
#' @examples
#' data(GrayFox_Hour)
#' lorelogram(GrayFox_Hour)
#'
#' @importFrom magrittr %>%
#' @export
lorelogram <- function(data, max_lag = 30, write_csv = TRUE, outDir = "", plot_LOR = TRUE, plot_title = "", plot_x_title = "Lag") {

  wd0 <- getwd()
  on.exit(setwd(wd0))

  if (class(data)!="data.frame") {
    stop("Data should be a data.frame",
         call. = FALSE)
  }

  # Remove rows (=cameras) with no detection and prepare data
  y <- data[rowSums(data[,2:ncol(data)], na.rm = TRUE) > 0,]
  y <- droplevels(y)
  y <- dplyr::rename(y, id=names(y[1]))


  # Determine all combinations of (current time, future times) for a sampling site up to max_lag.
  #
  # - V1 of x_cmb determines time point 1
  # - V2 of x_cmb determines time point 2
  #+ all_combinations
  max_delta_time <- max_lag

  if (ncol(y)>24*60*14) { #2-week data at minute time-interval
  # Create all combinations of minute-occasion up to 2 weeks of data
  x_cmb <- as.data.frame(arrangements::combinations(n=24*60*14, k=2, replace = FALSE))
  #head(x_cmb)

  # Now, calculate the time differences for all of these combinations and get rid of rows that include combinations of times where the difference in time <= max_delta_time.
  x_cmb <- x_cmb %>%
    dplyr::mutate(diff_time= x_cmb[,2] - x_cmb[,1]) %>% # column with time difference
    dplyr::filter(diff_time >= 0 & diff_time <= max_delta_time) %>% # filter combinations negative or too far apart
    dplyr::select(V1, V2)

  # Replicate to include more than 2 weeks of by-minute data
  x_cmb2 <- x_cmb
  for (i in 1:ceiling(ncol(y)/(24*60*14))){
    x_cmb_temp <- x_cmb + i*(24*60*14-max_delta_time)
    x_cmb2 <- rbind(x_cmb2,x_cmb_temp) %>% dplyr::filter(V1 <= ncol(y))
  }
  x_cmb <- x_cmb2
  rm(list=c("x_cmb2", "x_cmb_temp", "data"))

  # Remove duplicates and interval larger than time of the last occasion
  x_cmb <- x_cmb %>% dplyr::distinct() %>% dplyr::filter(V1 <= ncol(y) & V2 < ncol(y)) %>% dplyr::mutate(diff_time = V2 - V1)
  } else {
    x_cmb <- as.data.frame(arrangements::combinations(n=ncol(y), k=2, replace = FALSE))
    #head(x_cmb)

    # Now, calculate the time differences for all of these combinations and get rid of rows that include combinations of times where the difference in time <= max_delta_time.
    x_cmb <- x_cmb %>%
      dplyr::mutate(diff_time= x_cmb[,2] - x_cmb[,1]) %>% # column with time difference
      dplyr::filter(diff_time >= 0 & diff_time <= max_delta_time) %>% # filter combinations negative or too far apart
      dplyr::select(V1, V2)

    # Remove duplicates and interval larger than time last occasion
    x_cmb <- x_cmb %>% dplyr::distinct() %>% dplyr::filter(V1 <= ncol(y) & V2 < ncol(y)) %>% dplyr::mutate(diff_time = V2 - V1)
    }


  # #### Organize data: from wide format to long format
  #+ organize1
  # Organize data
  y2 <- tidyr::gather(y, time, value, -id)
  y2 <- dplyr::mutate(y2, time=as.numeric(substr(time,2,10)))
  #y3 <- y2 %>% filter(is.na(value)!=TRUE) #filter has to done after nest function
  y3 <- y2
  names(y3)[c(1,3)]<-c("id","y")
  y3$y <- as.numeric(as.character(y3$y))
  #head(y3)
  rm(list=c("y2"))

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
    Z<-Z %>% dplyr::mutate(time_diff=time2-time1)
    #head(Z)
    Z <- Z %>% dplyr::filter(!is.na(y1) & !is.na(y2))
    freq <- Z %>%
      dplyr::group_by(time_diff, y1, y2) %>%
      dplyr::summarise(count=n()) %>%
      tidyr::unite(col = y1y2, y1, y2, sep="", remove=FALSE) %>%
      dplyr::select(-y2) %>%  #remove info of second obs in time
      tidyr::spread(key=y1y2, value=count, fill=0)
    rm(list=c("y.t1","y.t2","Z"))
    freq
  }

  b <- lapply(1:n, myfunct)
  freq <- as.data.frame(data.table::rbindlist(b, fill = TRUE)) # from list to data.frame


# #### Organize data: Compile counts of pairwise 11, 01, 01, 00 for each time interval
#+ organize3
freq_tot <- freq %>%
  dplyr::group_by(time_diff, y1) %>%
  dplyr::summarise_all(dplyr::funs(sum))


# #### Calculate log odds ratios (direct calculation)
#+ Log_calculation
ORs<-freq_tot %>% dplyr::group_by(time_diff) %>%
  dplyr::summarize(or = sum(`11`)*sum(`00`)/(sum(`10`)*sum(`01`)),
                   se = sqrt(1/sum(`11`)+1/sum(`00`)+1/sum(`10`)+1/sum(`01`)))
ORs_all_minute <- ORs %>% dplyr::mutate(LORs=log(or),
                                        U_95_CI = log(or)+1.96*se,
                                        L_95_CI = log(or)-1.96*se)
ORs_all_minute[!is.na(ORs_all_minute$or) & ORs_all_minute$or==0,4:6] <- NA #replace values when lor values were undefined
LORs <- ORs_all_minute %>% dplyr::rename(., Lag=time_diff) %>% dplyr::select(c(Lag, LORs, U_95_CI, L_95_CI))

if (write_csv == TRUE) { # save results

  if(outDir != "" & !file.exists(outDir)){stop("outDir does not exist")}

  if (outDir == "") {
    filename <- paste("Log_odds_ratio_MaxLag_", as.character(max_delta_time),".csv", sep="")
    } else {
    filename <- paste(as.character(outDir), "/Log_odds_ratio_MaxLag_", as.character(max_delta_time),".csv", sep="")
    }
  write.csv(x = ORs_all_minute, file = filename)
}

if (plot_LOR == TRUE) { #plot lor values
plot_LOR_all <-  ggplot2::ggplot(ORs_all_minute, ggplot2::aes(x = time_diff)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = L_95_CI, ymax = U_95_CI, fill = "95 % CI"), alpha = 0.65) +
  ggplot2::geom_line(ggplot2::aes(y = LORs, color = "LORs"), size=0.25, linetype=1) + #
  ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype="solid")+
  ggplot2::scale_colour_manual("",values="black") +
  ggplot2::scale_fill_manual("",values="#4F6D7A") + #4F6D7A#F6511D
  ggplot2::labs(x = plot_x_title, y = "Log Odds Ratio", title = plot_title)+
  #coord_cartesian(ylim = c(-1,10), xlim=c(0,60))+ #
  ggplot2::theme_minimal()+
  ggplot2::theme(legend.justification = c(1, 1), legend.position = "none",
        axis.line.y = ggplot2::element_line(colour = 'black', linetype = 'solid'),
        axis.ticks.y = ggplot2::element_line(colour = 'black', linetype = 'solid'),
        axis.text = ggplot2::element_text(size=8),
        axis.title = ggplot2::element_text(size=10,face="bold"),
        panel.grid.minor.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(colour = 'grey', linetype = 'solid', size=0.25),
        panel.grid.minor.x = ggplot2::element_line(colour = 'grey', linetype = 'dashed', size=0.25))+
  ggplot2::scale_x_continuous(breaks=seq(0,max(ORs_all_minute$time_diff),10), labels=seq(0,max(ORs_all_minute$time_diff),10))
plot_LOR_all

  if (outDir == "") {
    filename <- paste("Log_odds_ratio_MaxLag_", as.character(max_delta_time),".jpg", sep="")
      } else {
    filename <- paste(as.character(outDir), "/Log_odds_ratio_MaxLag_", as.character(max_delta_time),".jpg", sep="")
      }
  ggplot2::ggsave(file = filename, plot_LOR_all, units = "cm", width = 30, height = 15)
  }
LORs
}
