
#'
#' Sample mean estimation of the mean
#'
#' Sample mean estimation of the population mean using the given sequence of the vertices.
#'
#'@param chain a sequence of node indices that represent the referral process.
#'@param trait.vec is a vector of size number of nodes in the network with 0-1 values that for example indicate a specific trait. 
#'
#'@return mean(trait.vec[chain], na.rm=TRUE) 
#'
#'@export
sample_mean  <- function(chain, trait.vec){
    return( mean(trait.vec[chain], na.rm=TRUE) )
}

#'
#' Volz-Heckathorn estimation of the mean
#'
#' Volz-Heckathorn estimation of the population mean using the given sequence of the vertices.
#'
#'@param chain a sequence of node indices that represent the referral process.
#'@param trait.vec is a vector of size number of nodes in the network with 0-1 values that for example indicate a specific trait. 
#'@param inc.prob is a vector of size number of nodes in the network with the corresponding stationary probabilities.
#'@references
#' Volz, Erik and Heckathorn, Douglas D (2008); ``Probability based estimation theory for respondent driven sampling'', Journal of official statistics.
#'
#'@details
#' Let x_1, x_2, \\cdots, x_n be a sequence of samples from a finite population. Furthermore, let \\pi_1, \\cdots, \\pi_N denote the stationary distribution. Then \\mu_\{VH\} = 1/(\\sum_\{i=1\}^n 1/\\pi_\{x_i\}) \\sum_\{i=1\}^n x_i/\\pi_\{x_i\}.
#'
#'@return the Volz-Heckathorn estimation of the population mean.
#'
#'@export
VH_mean  <- function(chain, trait.vec, inc.prob){
    stopifnot(length(inc.prob)==length(trait.vec))
    #NOTE: removing vertices with na trait.
    chain <- chain[ -which(chain %in% which(is.na(trait.vec)) ) ]

    if(length(chain) == 0)
        return(NA)

    stopifnot( length(which(is.na(trait.vec[chain]))) == 0)

    i.p <- inc.prob[chain] 
    return( sum(trait.vec[chain]/i.p, na.rm=TRUE)/sum(1/i.p) )
}



#'
#' Horvitz-Thompson estimation of the mean
#'
#' Horvitz-Thompson estimation of the population mean using the given sequence of the vertices.
#'
#'@param chain a sequence of node indices that represent the referral process.
#'@param trait.vec is a vector of size number of nodes in the network with 0-1 values that for example indicate a specific trait. 
#'@param inc.prob is a vector of size number of nodes in the network with the corresponding stationary probabilities.
#'
#'@references
#' Horvitz, D. G.; Thompson, D. J. (1952) ``A generalization of sampling without replacement from a finite universe'', Journal of the American Statistical Association,
#'
#'@details
#' Let x_1, x_2, \\cdots, x_n be a sequence of samples from a finite population of size N. Furthermore, let \\pi_1, \\cdots, \\pi_N be the inclusion of probabilities. Then \\mu_\{HT\} = 1/N*\\sum_\{i=1\}^n x_i/\\pi_\{x_i\}.
#'
#'@return the Horvitz-Thompson estimation of the population mean.
#'
#'@export
HT_mean  <- function(chain, trait.vec, inc.prob){
    stopifnot(length(inc.prob)==length(trait.vec))

    #NOTE: removing vertices with na trait.
    chain <- chain[ -which(chain %in% which(is.na(trait.vec)) ) ]

    if(length(chain) == 0)
        return(NA)

    stopifnot( length(which(is.na(trait.vec[chain]))) == 0)

    i.p <- inc.prob[chain] 
    return( mean(trait.vec[chain]/i.p, na.rm=TRUE)/length(trait.vec) )
}


