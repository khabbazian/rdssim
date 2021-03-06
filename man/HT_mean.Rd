% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RDSEstimators.R
\name{HT_mean}
\alias{HT_mean}
\title{Horvitz-Thompson estimation of the mean}
\usage{
HT_mean(chain, trait.vec, inclusion.probability)
}
\arguments{
\item{chain}{a sequence of node indices that represent the referral process.}

\item{trait.vec}{is a vector of size number of nodes in the network with 0-1 values that for example indicate a specific trait.}

\item{inclusion.probability}{is a vector of size number of nodes in the network with the corresponding stationary probabilities.}
}
\value{
the Horvitz-Thompson estimation of the population mean.
}
\description{
Horvitz-Thompson estimation of the population mean using the given sequence of the vertices.
}
\details{
Let x_1, x_2, \\cdots, x_n be a sequence of samples from a finite population of size N. Furthermore, let \\pi_1, \\cdots, \\pi_N be the inclusion of probabilities. Then \\mu_{HT} = 1/N*\\sum_{i=1}^n x_i/\\pi_{x_i}.
}
\references{
Horvitz, D. G.; Thompson, D. J. (1952) ``A generalization of sampling without replacement from a finite universe'', Journal of the American Statistical Association,
}
