#' Mixture CDF function with user-defined parameters
#'
#' This function computes the CDF of a mixture distribution.
#' @param q A numeric value where the CDF is evaluated.
#' @param mean1 The mean of the first normal distribution
#' @param sd1 The standard deviation corresponding to the first normal distribution.
#' @param mean2 The mean of the second normal distribution.
#' @param sd2 The second corresponding standard deviation.
#' @param mixprob The probability of the first distribution occurring.
#' @return A numeric value of the CDF at the original `q`.
#' @export

pnormmix <- function(q, mean1, sd1, mean2, sd2, mixprob) {
  (mixprob) * pnorm(q, mean = mean1, sd = sd1) + (1 - mixprob) * pnorm(q, mean = mean2, sd = sd2)
}

#' This function computes the PDF of a mixture distribution.
#' @param x a numeric value where the PDF is evaluated.
#' @param mean1 The mean of the first normal distribution
#' @param sd1 The standard deviation corresponding to the first normal distribution.
#' @param mean2 The mean of the second normal distribution.
#' @param sd2 The second corresponding standard deviation.
#' @param mixprob The probability of the first distribution occurring.
#' @return A numeric value for the PDF at the original `x`.
#' @export

dnormmix <- function(x, mean1, sd1, mean2, sd2, mixprob) {
  (mixprob) * dnorm(x, mean = mean1, sd = sd1) + (1 - mixprob) * dnorm(x, mean = mean2, sd = sd2)
}

#' This function generates a random sample from a mixture of two normal distributions.
#' @param n the number of desired samples.
#' @param mean1 The mean of the first normal distribution
#' @param sd1 The standard deviation corresponding to the first normal distribution.
#' @param mean2 The mean of the second normal distribution.
#' @param sd2 The second corresponding standard deviation.
#' @param mixprob The probability of the first distribution occurring.
#' @return a vector of the n samples.
#' @export


rnormmix <- function(n, mean1, sd1, mean2, sd2, mixprob) {
  coin <- sample(c("H", "T"), size = n, prob = c(mixprob, 1 - mixprob), replace = TRUE)
  X <- ifelse(coin == "H", rnorm(n, mean = mean1, sd = sd1), rnorm(n, mean = mean2, sd = sd2))
  return(X)
}

#' This function generates a random sample from a mixture of two normal distributions.
#' @param p a value for the quantile to be calculated.
#' @param mean1 The mean of the first normal distribution
#' @param sd1 The standard deviation corresponding to the first normal distribution.
#' @param mean2 The mean of the second normal distribution.
#' @param sd2 The second corresponding standard deviation.
#' @param mixprob The probability of the first distribution occurring.
#' @return the corresponding quantile to the input p
#' @export

qnormmix <- function(p, mean1, sd1, mean2, sd2, mixprob) {
  uniroot(function(x) pnormmix(x, mean1, sd1, mean2, sd2, mixprob) - p, interval = c(-10, 10), extendInt = "yes")$root
}

#' Vectorized quantile function
#'
#' This function is a vectorized version of the qnormmix function
#' @param p A vector of probabilities used to calculate different quantiles.
#' @param mean1 The mean of the first normal distribution
#' @param sd1 The standard deviation corresponding to the first normal distribution.
#' @param mean2 The mean of the second normal distribution.
#' @param sd2 The second corresponding standard deviation.
#' @param mixprob The probability of the first distribution occurring.
#' @return A vector of quantiles that corresponds to their respective probabilities of p
#' @export

myqnormmix <- Vectorize(qnormmix, vectorize.args = "p")



