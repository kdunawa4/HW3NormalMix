% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalmix.R
\name{rnormmix}
\alias{rnormmix}
\title{This function generates a random sample from a mixture of two normal distributions.}
\usage{
rnormmix(n, mean1 = 1, sd1 = 2, mean2 = 5, sd2 = 0.5, mixprob = 2/3)
}
\arguments{
\item{n}{the number of desired samples.}

\item{mean1}{Mean of the first normal distribution.}

\item{sd1}{Standard deviation of the first normal distribution.}

\item{mean2}{Mean of the second normal distribution.}

\item{sd2}{Standard deviation of the second normal distribution.}

\item{mixprob}{The probability of the first distribution occurring.}
}
\value{
a vector of the n samples.
}
\description{
This function generates a random sample from a mixture of two normal distributions.
}
