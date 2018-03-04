#'
#'Converts an adjacency matrix into the adjacency list
#'
#'@param A a symmetric unweighted adjacency matrix.
#'
#'@param ac.type indicates the anti clustering referral distribution types.
#'Type "Two" is "refer friends based on the number of friends they have that you
#'don't know", type "One" is "refer friends based on the number of friends you
#'have that don't know them" and type "Both" means type "One" and "Two".
#'
#'@export 
adj2list <- function(A, ac.type=c("Both", "Two", "One")){
    stopifnot( isSymmetric(A) )
    .Call('rdssim_adj2list', PACKAGE = 'rdssim', A, ac.type)
}



#'
#'Simulates referral chains/trees for Respondent-driven sampling
#'
#'@param adjList is the object \code{adj2list} returns. It has the adjacency
#'list of the input social network. It assumes that adjList represents an
#'undirected network. Refer to \code{adj2list}.
#'@param referral.type referral type. It can be "sRW", simple random walk, or
#'"acRW", anti-cluster random walk.
#'@param wRreplacement if TRUE, then the referral chain/tree will be collected with replacement.
#'@param nSamples number of samples to be collected.
#'@param nReferrals number of referral per participant (or node) 
#'@param seedNode the starting seed participant (or node).
#'@param rseed random seed for random number generators inside the code. It is useful for reproducing results.
#'
#'@examples
#'
#'generate.network <- function(N,con){
#'    A <- matrix(0,ncol=N,nrow=N) # initialize with zeros
#'    A[upper.tri(A)] <- rbinom(N*(N-1)/2,1,con) # fill up the upper triangle
#'    A[lower.tri(A)] <- t(A)[lower.tri(A)] # fill up the lower triangle to obtain a symmetric matrix
#'    A
#'}
#'
#'ER <- generate.network(300,0.6)
#'
#'adjL <- adj2list(ER)
#'res  <- rdssim(adjL, referral.type="sRW", wRreplacement = T, 
#'                 nSamples=100, nReferrals=3, seedNode=10, rseed=1)
#'
#'@references
#'Khabbazian, Mohammad, Bret Hanlon, Zoe Russek, and Karl Rohe. "Novel sampling
#'design for respondent-driven sampling." Electronic Journal of Statistics 11,
#'no. 2 (2017): 4769-4812.'
#'
#'@export 
rdssim <- function(adjList, referral.type=c("sRW","acRW"), wRreplacement=TRUE, nSamples, nReferrals, seedNode, rseed) {
    .Call('rdssim_rdssim_cpp', PACKAGE = 'rdssim', adjList$AdjList, adjL$AcAdjList, referral.type, wRreplacement, nSamples, nReferrals, seedNode, rseed)
}

