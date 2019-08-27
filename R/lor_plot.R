#' Save and customize lorelogram plot
#'
#'Estimates of pairwise log-odds ratios and associated 95\% confidence intervals for each lag between 1 and \code{max_lag} calculated by the \code{\link{lorelogram}} function are plotted against the temporal or spatial lag (x-axis). The plotted lorelogram provides a graphical description of how correlation between outcomes  changes as we increase the distance (in space or time) between sampling occassions. This function provides arguments for the customization of several aspects of the plot.
#'
#' @param data output of \code{\link{lorelogram}} function.
#' @param save_LOR_plot logical. Create a .jpg file of the results (default: FALSE)?
#' @param outDir character. Conditional on \code{save_LOR_plot} = TRUE. In which directory should the .jpg file be saved?
#' @param colour character. Color to fill the confidence interval band (default: "#0C71C9")
#' @param linetype character. Linetype for the average log-odds ratio values (default: "solid"). It accepts all the linetype available for ggplot2::geom_line.
#' @param title character. Title of the plot (default: NULL).
#' @param x_axis_title character. Title of x-axis of the plot (default: "Lag").
#' @param y_axis_title character. Title of y-axis of the plot (default: "Log-odds ratio").
#' @param ylim numeric. Vector of two values \code{c(y_min, y_max)} defining range for the y axis (default: NULL).
#' @param x_break numeric. Unit-lag distance between primary breaks in the x_axis (default: 10).
#' @param alpha numeric. Set transparency of the band describing the confidence interval. Should be a value between 0 and 1 (default: 1).
#' @return The function returns a plot of the estimates of pairwise log-odds ratios and associated 95\% confidence intervals for each lag between 1 and \code{max_lag}.
#' @details \code{data} must be a data.frame object containing the numeric output of the \code{\link{lorelogram}} function. Visual aspects of the lorelogram plot such as color and transparency of the confidence interval band, linetype of the curve representing the average log-odds ratios estimates, title and x- and y-axis labels, can be customized and the resulting plot can be saved in a folder of user choice.
#'
#' @examples
#' # import data and estimate log-odds ratio
#' data(GrayFox_Hour)
#' lor <- lorelogram(GrayFox_Hour, max_lag = 72, plot_LOR = FALSE)
#'
#' # basic plot
#' lor_plot(lor)
#'
#' # customized plot
#' lor_plot(lor, colour = "red", alpha = 0.7, title = "My lorelogram", x_break = 24, x_axis_title = "Time Lag (Hour)")
#'
#' @importFrom magrittr %>%
#' @export
lor_plot <- function(data, save_LOR_plot = FALSE,
                     outDir = "",
                     colour = "#0C71C9",
                     linetype = "solid",
                     title = "",
                     x_axis_title = "Lag",
                     y_axis_title = "Log-odds ratio",
                     ylim = NULL,
                     x_break = 10,
                     alpha = 1) {

  wd0 <- getwd()
  on.exit(setwd(wd0))



  plot_LOR_all <-  ggplot2::ggplot(data, ggplot2::aes(x = Lag)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = L_95_CI, ymax = U_95_CI, fill = "95 % CI"), alpha = alpha) +
      ggplot2::geom_line(ggplot2::aes(y = LORs, color = "LORs"), size=0.25, linetype=linetype) + #
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype="solid")+
      ggplot2::scale_colour_manual("",values="black") +
      ggplot2::scale_fill_manual("",values=colour) +
      ggplot2::labs(x = x_axis_title, y = y_axis_title, title = title)+{
        if (!is.null(ylim)) ggplot2::coord_cartesian(ylim = ylim)
      } +
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
       ggplot2::scale_x_continuous(breaks=seq(0,max(data$Lag),x_break), labels=seq(0,max(data$Lag),x_break))
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
