#' Identify lag to approximate independence in binary data
#'
#' This function identifies spatial or temporal intervals required to approximate independence in binary data by determining when the first derivative of the log-odds ratios with respect to \eqn{\Delta}\emph{t} reaches zero. It allows for non-zero asymptotic values due to among-site variation or cyclic patterns.
#'
#'
#' @param data output of \code{\link{lorelogram}} function.
#' @param n_knots numeric. Number of knots in the \code{\link[splines]{bs}}.
#' @param plot_results logical. Create a .jpg plot of the results (default: TRUE)?
#' @param outDir character. Directory into which .csv and plot file are saved.
#' @param plot_title character. Title of the plot (default: NULL).
#' @param plot_x_title character. Title of x-axis of the plot (default: Lag).
#' @return The function returns the minimum interval length necessary to approximate independence in the data.
#' @details \code{data} should be a data.frame containing the output of the \code{\link{lorelogram}} function.
#'
#' First, a smoothing cubic spline is fitted to the log-odds ratios values estimated by the \code{\link{lorelogram}} function at the series of spatial or temporal intervals of increasing length. Then, the first derivative of the spline curve is calculated.
#'
#' We recommend a visual inspection of the plot produced. When log-odds ratios present cyclic patterns (e.g., due to daily activity patterns affecting camera trap data) at the spatial or temporal scale used in the analysis, the \code{\link{lor_lag_to_indep}} can not be used to identify the minimum interval length that approximate independence in the data.
#'
#' @examples
#' data(GrayFox_Hour)
#' # cyclic pattern
#' lor <- lorelogram(GrayFox_Hour, max_lag = 120)
#' lor_lag_to_indep(lor)
#'
#' # non-cyclic pattern
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

a <- ggplot2::ggplot(lor, aes(x = Lag, y = LORs))+
  ggplot2::geom_point()+
  ggplot2::geom_line(aes(y = Y, color = "Y"), size = 1)+
  labs(x = "Lag", y = "Log-odds ratio", title = "Original fit")+
  scale_color_discrete(name = "", labels = "Cubic spline")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = '.'))+
  theme_classic()+
  theme(legend.position = c(0.8,0.95),
        legend.background = element_blank(),
        plot.title=element_text(hjust = 0.5, vjust = 0.5))
b <- ggplot2::ggplot(data = NULL, aes(x = dX, y = dY))+
  ggplot2::geom_line(size = 1)+
  labs(title = "Derivative")+
  theme_classic()+
  theme(plot.title=element_text(hjust = 0.5, vjust = 0.5))
c <- gridExtra::grid.arrange(a,b)

# Identify time-to-independence
which(abs(dY-0)==min(abs(dY-0))) #values closest to the zero

}
