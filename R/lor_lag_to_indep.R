#' Identify lag at which serial dependence is no longer present
#'
#' This function identifies the spatial or temporal lag at which serial dependence is no longer present in binary data by determining when the first derivative of the pairwise log-odds ratios with respect to \eqn{\Delta}\emph{t} reaches zero. It allows for non-zero asymptotic log-odds ratios due to non-stationarity in the mean process or among-site variation.
#'
#'
#' @param data output of \code{\link{lorelogram}} function.
#' @param n_knots numeric. Number of knots in the smoothing cubic spline (see splines::bs for details).
#' @param plot_results logical. Create a .jpg plot of the results (default: TRUE)?
#' @param outDir character. Directory into which .csv and plot file are saved.
#' @param plot_title character. Title of the plot (default: NULL).
#' @param plot_x_title character. Title of x-axis of the plot (default: Lag).
#' @return The function returns the minimum interval length necessary to approximate independence in the data.
#' @details \code{data} should be a data.frame containing the output of the \code{\link{lorelogram}} function.
#'
#' First, a cubic spline is fitted to the pairwise log-odds ratios estimated by the \code{\link{lorelogram}} function at the series of sampling occassions separated in time or space. Then, the first derivative of the spline curve is calculated. The lag at which the first derivative is closest to 0 (in absolute value) is returned. Users should visually inspect the plot of the first derivative function to ensure the lorelogram asymptotes and that this lag results in a minimum.
#'
#'
#' @examples
#' data(GrayFox_Hour)
#' lor <- lorelogram(GrayFox_Hour, max_lag = 30)
#' lor_lag_to_indep(lor)
#'
#'
#' @export
lor_lag_to_indep <- function(data, n_knots = 5, plot_results = TRUE, outDir = "", plot_title = "", plot_x_title = "Lag"){

lor <- data

# Fit cubic spline
# High number of knots to better approximate the LORs values
fit <- lm(lor$LORs~splines::bs(lor$Lag, knots = n_knots)) #bs used to fit a cubic spline.

# Calculate derivative
X <- data.frame(Lag=seq(min(lor$Lag),max(lor$Lag),by=1) ) # make an ordered sequence
Y <- predict(fit,newdata=X) # calculate predictions for that sequence
dY <- diff(Y)/diff(X$Lag)  # the derivative of the function
dX <- rowMeans(embed(X$Lag,2)) # centers the X values for plotting

if (plot_results == TRUE){
a <- ggplot2::ggplot(lor, ggplot2::aes(x = Lag, y = LORs))+
  ggplot2::geom_point()+
  ggplot2::geom_line(ggplot2::aes(y = Y, color = "Y"), size = 1)+
  ggplot2::labs(x = "Lag", y = "Log-odds ratio", title = "Original fit")+
  ggplot2::scale_color_discrete(name = "", labels = "Cubic spline")+
  ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = '.'))+
  ggplot2::theme_classic()+
  ggplot2::theme(legend.position = c(0.8,0.95),
        legend.background = ggplot2::element_blank(),
        plot.title=ggplot2::element_text(hjust = 0.5, vjust = 0.5))
b <- ggplot2::ggplot(data = NULL, ggplot2::aes(x = dX, y = dY))+
  ggplot2::geom_line(size = 1)+
  ggplot2::labs(title = "Derivative")+
  ggplot2::theme_classic()+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust = 0.5, vjust = 0.5))
c <- gridExtra::grid.arrange(a,b)
print(c)
}
# Identify time-to-independence
which(abs(dY-0)==min(abs(dY-0))) #values closest to zero

}
