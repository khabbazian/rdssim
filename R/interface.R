#'
#'Converts an adjacency matrix into the adjacency list
#'
#'@param A a symmetric unweighted adjacency matrix.
#'@export 
adj2list <- function(A) {
    stopifnot( isSymmetric(A) )
    .Call('rdssim_adj2list', PACKAGE = 'rdssim', A)
}

#'
#'Simulates referral chains/trees for Respondent-driven sampling
#'
#'@param adjlist adjacency list of the input social network. It assumes that adjlist represents an undirected network. Refer to \code{adj2list}.
#'@param rType referral type.
#'@param Markovian if TRUE and nReferrals > 1, then the sampling process follows Markov tree assumptions.
#'@param nSamples number of samples to be collected.
#'@param nReferrals number of referral per vertex.
#'@param seedNode the starting vertex.
#'@param rseed random seed for random number generators.
#'@export 
rdssim <- function(adjlist, rType, Markovian, nSamples, nReferrals = 1L, seedNode = 1L, rseed = 1L) {
    .Call('rdssim_rdssim_cpp', PACKAGE = 'rdssim', adjlist, rType, Markovian, nSamples, nReferrals, seedNode, rseed)
}

