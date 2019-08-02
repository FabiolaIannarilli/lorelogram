#' Plot pairwise log-odds ratios at intervals of increasing length
#'
#' The estimates of pairwise log-odds ratios and associated 95\% confidence intervals for each lag between 1 and \code{max_lag} calculated by the \code{\link{lorelogram}} function are plotted against (spatial or temporal) intervals of increasing length. The plot of the estimates versus lags provides a graphical description of how correlation between outcomes x-lag apart changes at the increase of the distance (in space or time) between the sampling replicates and it allows to identify which spatial or temporal interval is necessary to approximate independence in binary data.
#'
#' @param data output of \code{\link{lorelogram}} function.
#' @param save_LOR_plot logical. Create a .jpg plot of the results (default: FALSE)?
#' @param outDir character. Directory into which .csv and plot file are saved.
#' @param colour character. Color to fill the confidence interval band (defaul: "#0C71C9")
#' @param plot_title character. Title of the plot (default: NULL).
#' @param plot_x_title character. Title of x-axis of the plot (default: Lag).
#' @return The fucntion returns a plot of the estimates of pairwise log-odds ratios and associated 95\% confidence intervals for each lag between 1 and \code{max_lag}.
#' @details \code{data} should resemble a binary detection/nondetection history matrix such that provided as an output by the function \code{\link[camtrapR:detectionHistory]{camtrapR::detectionHistory}}. \code{\link{lorelogram}} can handle NAs in \code{data}.
#'
#' @examples
#' data(GrayFox_Hour)
#' lor <- lorelogram(GrayFox_Hour, max_lag = 120)
#' lor_plot(lor)
#'
#' @importFrom magrittr %>%
#' @export
lor_plot <- function(data, save_LOR_plot = FALSE, outDir = "", colour = "#0C71C9", plot_title = "", plot_x_title = "Lag") {

  wd0 <- getwd()
  on.exit(setwd(wd0))



  plot_LOR_all <-  ggplot2::ggplot(data, ggplot2::aes(x = Lag)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = L_95_CI, ymax = U_95_CI, fill = "95 % CI"), alpha = 0.65) +
      ggplot2::geom_line(ggplot2::aes(y = LORs, color = "LORs"), size=0.25, linetype=1) + #
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype="solid")+
      ggplot2::scale_colour_manual("",values="black") +
      ggplot2::scale_fill_manual("",values="#0C71C9") +
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
      ggplot2::scale_x_continuous(breaks=seq(0,max(data$Lag),10), labels=seq(0,max(data$Lag),10))
    plot_LOR_all

  if (save_LOR_plot == TRUE) { #save plot
    if (outDir == "") {
      filename <- paste("Log_odds_ratio_MaxLag_", as.character(max(data$Lag)),".jpg", sep="")
    } else {
      filename <- paste(as.character(outDir), "/Log_odds_ratio_MaxLag_", as.character(max(data$Lag)),".jpg", sep="")
    }
    ggplot2::ggsave(file = filename, plot_LOR_all, units = "cm", width = 30, height = 15)
  }
plot_LOR_all
}
