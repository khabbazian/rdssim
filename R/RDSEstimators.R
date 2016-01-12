

#'@param chain a sequence of node indices that represent the referral process.
#'@param trait.vec is a vector of size number of nodes in the network with corresponding trait values.
sample_mean  <- function(chain, trait.vec){
    return( mean(trait.vec[chain[,2]], na.rm=TRUE) )
}



#'@param chain a sequence of node indices that represent the referral process.
#'@param trait.vec is a vector of size number of nodes in the network with corresponding trait values.
#'@param inclusion.probability is a vector of size number of nodes in the network with corresponding stationary probabilities.
HT_mean  <- function(chain, trait.vec, inclusion.probability){
    stopifnot(length(inclusion.probability)==length(trait.vec))
    #NOTE: just the node sequence
    chain <- chain[,2]
    #NOTE: removing vertices with na trait.
    chain <- chain[ -which(chain %in% which(is.na(trait.vec)) ) ]

    if(length(chain) == 0)
        return(NA)

    stopifnot( length(which(is.na(trait.vec[chain]))) == 0)

    i.p <- inclusion.probability[chain] 
    return( sum(trait.vec[chain]/i.p, na.rm=TRUE)/sum(1/i.p) )
}



#'@param chain a sequence of node indices that represent the referral process.
#'@param trait.vec is a vector of size number of nodes in the network with corresponding trait values.
#'@param inclusion.probability is a vector of size number of nodes in the network with corresponding stationary probabilities.
HH_mean  <- function(chain, trait.vec, inclusion.probability){
    stopifnot(length(inclusion.probability)==length(trait.vec))
    #NOTE: just the node sequence
    chain <- chain[,2]
    #NOTE: removing vertices with na trait.
    chain <- chain[ -which(chain %in% which(is.na(trait.vec)) ) ]

    if(length(chain) == 0)
        return(NA)

    stopifnot( length(which(is.na(trait.vec[chain]))) == 0)

    i.p <- inclusion.probability[chain] 
    return( mean(trait.vec[chain]/i.p, na.rm=TRUE)/length(trait.vec) )
}


