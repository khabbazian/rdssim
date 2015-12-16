

#' Sample mean estimator
#'
#'@param chain a sequence of node indices that represent the referral process.
#'@param trait.vec is a vector of size number of nodes in the network with corresponding trait values.
#'
#'@export
#'
sample_mean  <- function(chain, trait.vec){
    return( mean(trait.vec[chain[,2]], na.rm=TRUE) )
}



#' HT estimator
#'
#'@param chain a sequence of node indices that represent the referral process.
#'@param trait.vec is a vector of size number of nodes in the network with corresponding trait values.
#'@param inclusion.probability is a vector of size number of nodes in the network with corresponding stationary probabilities.
#'
#'@export
#'
HT_mean  <- function(chain, trait.vec, inclusion.probability){
    stopifnot(length(inclusion.probability)==length(trait.vec))
    i.p <- inclusion.probability[chain[,2]] 
    return( sum(trait.vec[chain[,2]]/i.p, na.rm=TRUE)/sum(1/i.p) )
}



#' HH estimator
#'
#'@param chain a sequence of node indices that represent the referral process.
#'@param trait.vec is a vector of size number of nodes in the network with corresponding trait values.
#'@param inclusion.probability is a vector of size number of nodes in the network with corresponding stationary probabilities.
#'
#'@export
#'
HH_mean  <- function(chain, trait.vec, inclusion.probability){
    stopifnot(length(inclusion.probability)==length(trait.vec))
    chain <- chain[,2]
    i.p <- inclusion.probability[chain] 
    return( mean(trait.vec[chain]/i.p, na.rm=TRUE)/length(trait.vec) )
}
