#' One-Sample t-Test with Wallace Critical Values and Hastings P-Values
#'
#' Performs a comprehensive one-sample t-test using Wallace approximation for critical values 
#' and 4-term Hastings approximation for p-values. Includes normality testing, outlier detection,
#' and diagnostic plots (histogram, Q-Q plot, boxplot).
#'
#' @param x Character string. Name of the numeric column in `data` to test.
#' @param mu Numeric. Hypothesized population mean under the null hypothesis (H0: mu = mu0).
#' @param ha Character string. Alternative hypothesis. One of "less" (H1: mu < mu0), 
#'   "greater" (H1: mu > mu0), or "not equal" (H1: mu != mu0).
#' @param alpha Numeric. Significance level (0 < alpha < 1). Must be explicitly provided.
#' @param data Data frame containing the column specified in `x`.
#' @param plot Logical. If TRUE (default), display diagnostic plots: Q-Q plot, histogram 
#'   with normal curve, and boxplot.
#'
#' @return List of class "oneMeanTTest" containing these fields:
#' \itemize{
#'   \item variable: Tested variable name
#'   \item n: Sample size (after removing NAs)
#'   \item sample_mean: Sample mean
#'   \item sample_sd: Sample standard deviation  
#'   \item t_statistic: t-test statistic
#'   \item p_value: Hastings p-value
#'   \item confidence_interval: Confidence bounds
#'   \item decision: "REJECT" or "FAIL TO REJECT"
#'   \item is_normal: Normality assumption met?
#'   \item outliers: Detected outliers
#' }
#'
#' @details 
#' Critical values use Wallace approximation: \deqn{t_\alpha = \sqrt{\nu (\exp(c \cdot z_\alpha^2 / \nu) - 1)}}
#' where \eqn{c = ((8\nu+3)/(8\nu+1))^2} and \eqn{\nu = n-1}.
#'
#' P-values use 4-term Hastings approximation for computational efficiency.
#'
#' @section Assumptions:
#' 1. Normality (Shapiro-Wilk test, p > alpha)
#' 2. No outliers (1.5 * IQR rule)
#'
#' @examples
#' # Example 1: Earnings data (suspected to be below $75k)
#' set.seed(181)
#' earnings_data <- data.frame(
#'   Earnings = round(rnorm(50, mean = 65, sd = 20), 0)
#' )
#' result <- oneMeanTTest("Earnings", mu = 75, ha = "less", 
#'                        alpha = 0.05, data = earnings_data)
#' print(result)
#'
#' # Example 2: Test scores (check if above 70)
#' set.seed(456)
#' scores <- data.frame(
#'   Score = round(runif(40, min = 50, max = 95), 1)
#' )
#' oneMeanTTest("Score", mu = 70, ha = "greater", 
#'              alpha = 0.05, data = scores, plot = FALSE)
#'
#' # Example 3: Temperature readings (two-tailed test)
#' temp_data <- data.frame(
#'   Temperature = c(98.2, 98.6, 97.8, 99.1, 98.4, 98.7, 97.9, 
#'                   98.5, 98.3, 98.8, 98.1, 98.6)
#' )
#' oneMeanTTest("Temperature", mu = 98.6, ha = "not equal", 
#'              alpha = 0.05, data = temp_data)
#'
#' @export
oneMeanTTest <- function(x, mu, ha, alpha, data, plot = TRUE) {
  
  # Input validation
  if (is.null(data)) stop("Error: Data cannot be null.")
  if (!is.data.frame(data)) stop("Error: Data must be a data frame.")
  
  if (!is.character(x) || length(x) != 1) stop("Error: x must be a single column name (character).")
  if (!(x %in% names(data))) stop("Error: x is not a column in data.")
  
  if (!is.numeric(mu)) stop("Error: Hypothesized mean must be numeric.")
  if(missing(alpha)){
    stop("Error: 'alpha' or significance level must be provided")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("Error: Alpha must be between 0 and 1.")
  
  valid_alternatives <- c("less", "greater", "not equal")
  if (!(ha %in% valid_alternatives)) {
    stop("Error: ha must be 'less', 'greater', or 'not equal'.")
  }
  
  # Extract and validate data
  col_data <- data[[x]]
  
  if (all(is.na(col_data))) stop("Error: All values are NA.")
  col_data <- na.omit(col_data)
  
  if (!is.numeric(col_data)) stop("Error: Column must contain numeric values.")
  if (length(col_data) < 2) stop("Error: Need at least two observations for t-test.")
  if (sd(col_data) == 0) stop("Error: All observations identical; t-test cannot be performed.")
  
  # Assumption checks: Normality (Shapiro-Wilk)
  shapiro_test <- shapiro.test(col_data)
  data_is_normal <- shapiro_test$p.value > alpha
  
  # Assumption checks: Outliers (1.5 * IQR Rule) 
  Q1 <- quantile(col_data, 0.25)
  Q3 <- quantile(col_data, 0.75)
  IQR_val <- Q3 - Q1
  
  lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val
  
  outliers <- col_data[col_data < lower_bound | col_data > upper_bound]
  has_outliers <- length(outliers) > 0
  
  # Diagnostic plots
  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    layout(matrix(c(1, 2, 3), nrow = 1))
    
    # Q-Q Plot (Normality check)
    qqnorm(col_data, main = "Q-Q Plot", col = ifelse(data_is_normal, "orange", "red"))
    qqline(col_data, col = "darkblue", lwd = 2)
    mtext(paste("Shapiro-Wilk p =", round(shapiro_test$p.value, 4)))
    
    # Histogram with normal curve overlay
    hist_output <- hist(col_data, col = "darkblue", main = "Histogram with Normal Curve",
                        xlab = "Data Value", ylab = "Frequency", border = "white")
    
    m <- mean(col_data)
    s <- sd(col_data)
    
    # Normal curve (scaled to match histogram counts)
    x_vals <- seq(min(col_data), max(col_data), length = 200)
    y_vals <- dnorm(x_vals, mean = m, sd = s)
    y_vals <- y_vals * diff(hist_output$mids[1:2]) * length(col_data)
    
    lines(x_vals, y_vals, col = "red", lwd = 2)
    
    legend("topright", legend = c("Normal curve", paste("Normality:", ifelse(data_is_normal, "PASS", "FAIL"))),
           col = c("red", NA), lty = c(1, NA), lwd = c(2, NA), bty = "n")
    
    # Boxplot (Outlier check)
    boxplot(col_data, main = "Outlier Check", col = ifelse(has_outliers, "red", "darkblue"))
    if(has_outliers) {
      text(x = 1, y = outliers, labels = round(outliers, 2), pos = 4, col = "red", cex = 0.8)
    }
  }
  
  # Compute test statistic
  n <- length(col_data)
  xbar <- mean(col_data)
  sd_val <- sd(col_data)
  t_stat <- (xbar - mu) / (sd_val / sqrt(n))
  df <- n - 1
  
  # Wallace critical value approximation (Z to T conversion)
  z_input <- switch(ha,
                    "greater"   = qnorm(1 - alpha),
                    "less"      = qnorm(alpha),
                    "not equal" = qnorm(1 - alpha/2))
  
  term_crit <- (z_input^2 / df) * ((8*df + 3) / (8*df + 1))^2
  t_critical <- sqrt(df * (exp(term_crit) - 1))
  
  if(ha == "less") t_critical <- -t_critical
  
  # Hastings p-value approximation (4-term)
  term_p <- df * log(1 + (t_stat^2 / df))
  correction_factor <- (8*df + 1) / (8*df + 3)
  z_stat_equiv <- sqrt(term_p) * correction_factor
  
  if (t_stat < 0) z_stat_equiv <- -z_stat_equiv
  
  # Hastings cumulative distribution function
  Phi_hastings <- function(z) {
    if (z < 0) return(1 - Phi_hastings(-z))
    a1 <- 0.196854
    a2 <- 0.115194
    a3 <- 0.000344
    a4 <- 0.019527
    
    t <- 1 + a1*z + a2*z^2 + a3*z^3 + a4*z^4
    return(1 - 0.5 * t^(-4))
  }
  
  # Compute p-value based on alternative hypothesis
  p_value <- switch(ha,
                    "not equal" = 2 * (1 - Phi_hastings(abs(t_stat))),
                    "less"      = Phi_hastings(t_stat),
                    "greater"   = 1 - Phi_hastings(t_stat))
  
  # Confidence interval construction
  se <- sd_val / sqrt(n)
  margin <- t_critical * se
  
  if (ha == "not equal") {
    ci <- c(lower = xbar - margin, upper = xbar + margin)
    ci_text <- paste("between", round(ci[1], 4), "and", round(ci[2], 4))
    
  } else if (ha == "greater") {
    ci <- c(lower = xbar - margin, upper = Inf)
    ci_text <- paste("greater than", round(ci[1], 4))
    
  } else { # less
    ci <- c(lower = -Inf, upper = xbar - margin)
    ci_text <- paste("less than", round(ci[2], 4))
  }
  
  # Decision and interpretation
  decision <- ifelse(p_value < alpha, "REJECT the null hypothesis (Ho)", "FAIL TO REJECT the null hypothesis (Ho)")
  
  conf_level <- (1 - alpha) * 100
  
  interpretation <- paste0(
    "At alpha = ", alpha, ", we ", decision, ". The sample mean (", round(xbar, 4), ") is ", 
    ifelse(p_value < alpha, "significantly ", "not significantly "), ha, 
    ifelse(ha == "not equal", " to ", " than "), mu, ".\n",
    "We are ", conf_level, "% confident that the true population mean is ", ci_text, "."
  )
  
  # Return S3 object
  result <- list(
    variable = x,
    n = n,
    sample_mean = xbar,
    sample_sd = sd_val,
    hypothesized_mean = mu,
    t_statistic = t_stat,
    df = df,
    critical_value = t_critical,
    p_value = p_value,
    confidence_interval = ci,
    alternative = ha,
    alpha = alpha,
    decision = decision,
    interpretation = interpretation,
    shapiro_stat = shapiro_test$statistic, 
    shapiro_p = shapiro_test$p.value, 
    is_normal = data_is_normal,
    outlier_count = length(outliers), 
    outliers = outliers
  )
  
  class(result) <- "oneMeanTTest"
  return(result)
}

#' Print Method for oneMeanTTest Objects
#'
#' Pretty printing for oneMeanTTest results with assumption checks and interpretation.
#'
#' @param x An object of class "oneMeanTTest".
#' @param ... Additional arguments (ignored).
#'
#' @keywords internal
print.oneMeanTTest <- function(x, ...) {
  cat("\n")
  cat("    One-Sample t-test\n")
  cat("------------------------------------------------------\n")
  cat("Variable:   ", x$variable, "\n")
  cat("n:          ", x$n, "\n")
  cat("Mean:       ", round(x$sample_mean, 4), "\n")
  cat("SD:         ", round(x$sample_sd, 4), "\n")
  cat("------------------------------------------------------\n")
  cat("Null Hyp:    mu =", x$hypothesized_mean, "\n")
  cat("Alt Hyp:     mu", switch(x$alternative, "less"="<", "greater"=">", "not equal"="!="), 
      x$hypothesized_mean, "\n")
  cat("------------------------------------------------------\n")
  cat("ASSUMPTION CHECKS:\n")
  
  # Normality check
  cat("1. Normality (Shapiro-Wilk): ")
  if(x$is_normal) cat("SATISFIED (p =", format.pval(x$shapiro_p, digits=3), ")\n")
  else cat("NOT SATISFIED (p =", format.pval(x$shapiro_p, digits=3), ") - Consider Wilcoxon Test\n")
  
  # Outlier check
  cat("2. Outliers (1.5 IQR Rule):  ")
  if(x$outlier_count == 0) {
    cat("SATISFIED (0 outliers detected)\n")
  } else {
    cat("WARNING (", x$outlier_count, " detected: ", paste(round(x$outliers, 2), collapse=", "), ")\n", sep="")
  }
  
  cat("------------------------------------------------------\n")
  cat("t-stat:     ", round(x$t_statistic, 4), "\n")
  cat("df:         ", x$df, "\n")
  cat("Crit Val:   ", round(x$critical_value, 4), "\n")
  cat("P-value:    ", format.pval(x$p_value, digits=4), "\n")
  cat("------------------------------------------------------\n")
  
  # Display confidence interval
  ci_lower_str <- ifelse(is.infinite(x$confidence_interval[1]), "-Inf", round(x$confidence_interval[1], 4))
  ci_upper_str <- ifelse(is.infinite(x$confidence_interval[2]), "Inf", round(x$confidence_interval[2], 4))
  
  cat(paste0((1 - x$alpha)*100, "% CI:     "), 
      "[", ci_lower_str, ", ", ci_upper_str, "]\n")
  
  cat("Decision:   ", x$decision, "\n")
  cat("\nInterpretation:\n")
  cat(x$interpretation, "\n")
  cat("\n")
}
