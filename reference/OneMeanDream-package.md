# OneMeanDream: Enhanced One-Sample t-Test with Diagnostics

Performs one-sample t-tests using Wallace critical value approximation
and Hastings p-value approximation. Includes comprehensive assumption
checks (normality via Shapiro-Wilk, outlier detection) and diagnostic
plots.

## Main Functions

- [`oneMeanTTest`](https://kingssibal.github.io/OneMeanDream/reference/oneMeanTTest.md) -
  Perform one-sample t-test with diagnostics

- [`print.oneMeanTTest`](https://kingssibal.github.io/OneMeanDream/reference/print.oneMeanTTest.md) -
  Pretty print results

## Key Features

- Wallace approximation for critical values

- 4-term Hastings approximation for p-values

- Shapiro-Wilk normality test

- 1.5 \* IQR outlier detection

- Diagnostic plots (Q-Q plot, histogram, boxplot)

- Detailed interpretation of results

## Getting Started

Load the package and run a basic test:

[`library(OneMeanDream)`](https://kingssibal.github.io/OneMeanDream/)

`# Sample data` `set.seed(123)`
`data <- data.frame(x = rnorm(30, mean = 10, sd = 2))`

`# Run test`
`result <- oneMeanTTest("x", mu = 10, ha = "not equal", alpha = 0.05, data = data)`

`# View results` `print(result)`

## See also

Useful links:

- <https://kingssibal.github.io/OneMeanDream/>

## Author

**Maintainer**: Reijel Agub <rcagub@up.edu.ph>

Authors:

- Arianna Tagulinao

- Cindy Barcena

- Francesca Anne Baclig

- Josiah Tan

- Junel Rey

- Keano Pastrana

- Matt Andrei Miclat

- Mica Gonzales

- Zan Eros Policarpio
