% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RDSEstimators.R
\name{cmp_cov}
\alias{cmp_cov}
\title{Computes the covariance between the samples in the Markov chain/tree}
\usage{
cmp_cov(A, covariates, distance = 1)
}
\arguments{
\item{A}{adjacency matrix of the undirected network.}

\item{covariates}{the characteristic of interest for each node in the network.}

\item{distance}{distance between the random variables in the Markov chain/tree.}
}
\description{
Computes the covariance between the samples in the Markov chain/tree using
eigenvalues and eigenvectors of the Laplacian matrix of the social network.
}
\details{
Let \\pi, y, and d be the stationary distribution, the covariate and the
distance respectively. Then cov(y(X_i),y(X_{i+d})) = \\sum_{i=2}^N <y,f_i>^2_\\pi
\\lambda_i^d.
}
