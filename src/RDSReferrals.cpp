//
// Under GNU (>=2) license 
//
// flags for Rcpp 
// [[Rcpp::depends("RcppEigen")]]
// [[Rcpp::plugins(cpp11)]]

#include "types.h"

//NOTE: SIMP_RW:  A simple Markovian random walk.
//NOTE: AC_RW:  AC random walk choses uniformly from the set of nodes form wedge
//NOTE: with current node.


// returns a random integer X in the set {l, l+1, ..., u}.
int random_int(const int l, const int u){
    RASSERT( l <= u );
    double uR = l + (u+1 - l)*(double) rand()/((double) RAND_MAX + 1);
    return (int) uR;
}

enum ReferralType{SIMP_RW, AC_RW};

inline bool are_not_neighbor(const AdjList &adjlist, const int node1, const int node2){

    if(node1 == node2) 
        return 0;

    const auto l = adjlist[node1];
    return !(std::binary_search(l.begin(), l.end(), node2));
}

inline bool are_neighbor(const AdjList &adjlist, const int node1, const int node2){

    if(node1 == node2) 
        return 1;

    const auto l = adjlist[node1];
    return std::binary_search(l.begin(), l.end(), node2);
}

template<int nReferrals>
vector<Referral> refer_next_node(const AdjList &, const int, const int, const ReferralType);

template<int nReferrals>
vector<Referral> refer_next_node_ac_rw(const AdjList &adjlist, const int cNode, const int prevNode, 
        const ReferralType rt){

    double weight=1;
    const auto list = adjlist[cNode];

    NodePairVector possibleReferrals;
    possibleReferrals.reserve(1024);

    for(auto item1:list)
        for(auto item2:list)
            if ( are_not_neighbor(adjlist, item1, item2) )
                possibleReferrals.push_back( NodePair(item1, item2) );

    const int typeOneNCandids = possibleReferrals.size();

    if(nReferrals == 1){
        for(auto item1:list)
            for(auto item2:adjlist[item1])
                if ( are_not_neighbor(adjlist, item2, cNode) )
                    possibleReferrals.push_back( NodePair(item1, item2) );
    } else { // refers a node that is not connected to the prev node
        for(auto item1:list)
            if ( are_not_neighbor(adjlist, item1, prevNode) )
                possibleReferrals.push_back( NodePair(item1, prevNode) );
    }

    const int nCandids = possibleReferrals.size() ;
    const int typeTwoNCandids = nCandids - typeOneNCandids;

    if( !nCandids ){//Handling special cases
        Rf_warning(("cN "+to_string(cNode)+". Zero candid list!").c_str());
        return refer_next_node<nReferrals>(adjlist, cNode, prevNode, (ReferralType) 0);
    }

    if( nReferrals == 1 || nCandids == 1 ){ //single referral for the Markov chain/tree simulation

        const int id = random_int(0, nCandids-1);

        RASSERT(id < possibleReferrals.size() );
        const int nextNode = get<0>( possibleReferrals[id] );
        return vector<Referral> { make_tuple(nextNode, weight) };

    } else if( nReferrals == 2 ){

        if( typeOneNCandids > 0 && typeTwoNCandids > 0){

            const int idOne = random_int(0, typeOneNCandids-1);
            const int idTwo = typeOneNCandids + random_int(0, typeTwoNCandids-1);
            const auto node1 = get<0>( possibleReferrals[idOne] );
            const auto node2 = get<0>( possibleReferrals[idTwo] );
            if( node1 != node2 )
                return vector<Referral> {make_tuple(node1,weight), make_tuple(node2,weight)};
            else
                return vector<Referral> {make_tuple(node1,weight)};

        } else if( typeOneNCandids == 0 ){
            const int idTwo = typeOneNCandids + random_int(0, typeTwoNCandids-1);
            const auto node2 = get<0>( possibleReferrals[idTwo] );
            return vector<Referral> {make_tuple(node2,weight)};
        } else if( typeTwoNCandids == 0 ){
            const int idOne = random_int(0, typeOneNCandids-1);
            const auto node1 = get<0>( possibleReferrals[idOne] );
            const auto node2 = get<1>( possibleReferrals[idOne] );
            return vector<Referral> {make_tuple(node1,weight), make_tuple(node2,weight)};
        } else
            RASSERT(0);

    } else if( nReferrals == 3 ){

        vector<Referral> refVec;
        int node1, node2;
        if( typeOneNCandids > 0 ){
            const int idOne = random_int(0, typeOneNCandids-1);
            node1 = get<0>( possibleReferrals[idOne] );
            node2 = get<1>( possibleReferrals[idOne] );

            refVec.push_back( make_tuple(node1, weight) );
            if(node1 != node2 )
                refVec.push_back( make_tuple(node2, weight) );
        }

        if( typeTwoNCandids > 0 ){
            const int idTwo  = typeOneNCandids + random_int(0, typeTwoNCandids-1);
            const auto node3 = get<0>( possibleReferrals[idTwo] );

            if(node3 != node1 && node3 != node2)
                refVec.push_back( make_tuple(node3, weight) );
        }

        RASSERT(refVec.size() > 0);
        return refVec;
    } 

}


template<int nReferrals>
vector<Referral> refer_next_node(const AdjList &adjlist, const int cNode, const int prevNode, 
        const ReferralType rt){

    RASSERT( cNode < adjlist.size() && cNode >= 0);
    RASSERT( prevNode < adjlist.size() && prevNode >= 0);
    RASSERT( nReferrals < 4 && nReferrals > 0 );

    if( rt == AC_RW ) 
        return refer_next_node_ac_rw <nReferrals>(adjlist, cNode, prevNode, rt);


    const auto list      = adjlist[cNode];
    const int nNeighbors = list.size();
    const double weight  = 1/(double) nNeighbors;

    vector<Referral> refVec;

    const int firstIdx = random_int(0, nNeighbors-1);
    const int firstNode = list[ firstIdx ];

    refVec.push_back( make_tuple(firstNode, weight) );

    if (nReferrals == 1 || nNeighbors < 2) 
        return refVec;

    //NOTE: scenario for referring the second node; a w/o replacement sampling from neighbors.
    const int secondIdx  = random_int(0, nNeighbors-1);
    const int secondNode = list[ secondIdx ];

    if( secondNode != firstNode )
        refVec.push_back( make_tuple( secondNode, weight) );

    if (nReferrals == 2 || nNeighbors < 3) 
        return refVec;

    const int thirdIdx = random_int(0, nNeighbors-1);
    const int thirdNode = list[ thirdIdx ];

    if( thirdNode != firstNode && thirdNode != secondNode )
        refVec.push_back( make_tuple( thirdNode, weight) );

    return refVec;

} 

inline void fill_log(Matrix &logs, const int currentNode, const Referral nextReferral, const int wave, const int l){
        //NOTE: currentNode; It has to be modified to R array format that indexing starts from 1
        logs(l,0) = currentNode + 1;
        logs(l,1) = get<0>( nextReferral ) + 1;
        //NOTE: the probability of choosing the next node from the candid list
        logs(l,2) = get<1>( nextReferral );
        //NOTE: wave 
        logs(l,3) = wave;
}

Matrix sim_referral_chain(const AdjList &adjlist,
        const int chainLength, const int seedNode, const ReferralType rt, 
        const unsigned int rseed, Matrix &logs){

    RASSERT(logs.nrow() == chainLength);
    RASSERT(logs.ncol() == 4);

    if(adjlist[seedNode].size() == 0)
        Rf_error("The seed node is an isolated node!");

    int currentNode = seedNode, previousNode = seedNode;

    for(int level = 0; level<chainLength; ++level){

        const int nNeighbors = adjlist[currentNode].size();
        if (nNeighbors == 0)
            Rf_error("Zero neighbors! That should not happen.");

        const auto nextReferralVec = refer_next_node<1>(adjlist, currentNode, previousNode, rt);
        const auto nextReferral    = nextReferralVec[0];

        //NOTE: I assume that the node that enters the sampling is referred directly by a participant.
        RASSERT( are_neighbor(adjlist, currentNode, get<0>( nextReferral ) ) );

        fill_log(logs, currentNode, nextReferral, level, level); 

        previousNode = currentNode;
        currentNode  = get<0>(nextReferral);
    }
    return logs;
}


Matrix sim_referral_tree(const AdjList &adjlist,
        const int nSamples, const int nReferrals,
        const int seedNode, 
        const ReferralType rt, 
        const bool Markovian,
        const unsigned int rseed, Matrix &logs){

    if (adjlist[seedNode].size() == 0)
        Rf_error("The seed node is an isolated node!");

    RASSERT(logs.nrow() == nSamples);
    RASSERT(logs.ncol() == 4);

    int currentNode = seedNode, previousNode = seedNode, wave=0;
    double weight = 1;

    std::set<int> visitedNodes;

    //node, parent, weight, wave number
    typedef tuple<int, int, double, int> QObject;
    queue<QObject> toBeVisited;
    toBeVisited.push( QObject(currentNode, previousNode, weight, wave) );

    int counter=0;
    for(; counter<nSamples; ++counter){

        const int nNeighbors = adjlist[currentNode].size();
        RASSERT( nNeighbors > 0 );
        if( toBeVisited.size() == 0 )
            break;

        const auto qObj = toBeVisited.front(); 
        toBeVisited.pop();

        currentNode  = get<0>(qObj);
        previousNode = get<1>(qObj);
        weight       = get<2>(qObj);
        wave         = get<3>(qObj);


        if( visitedNodes.find(currentNode) != visitedNodes.end() ){
            --counter;
            continue;
        }
        visitedNodes.insert(currentNode);

        fill_log(logs, previousNode, make_tuple(currentNode, weight), wave, counter); 

        if( Markovian == true ){
            for(int i=0; i < nReferrals; ++i){
                const auto nextReferralVec = refer_next_node<1>(adjlist, currentNode, previousNode, rt);
                const auto nextReferral    = nextReferralVec[0];
                toBeVisited.push( QObject( get<0>(nextReferral), currentNode, get<1>( nextReferral ), wave+1) );
            }
        } else {
            
            vector<Referral>  nextReferralVec;
            if(nReferrals == 1)
                nextReferralVec = refer_next_node<1>(adjlist, currentNode, previousNode, rt);
            else if(nReferrals == 2)
                nextReferralVec = refer_next_node<2>(adjlist, currentNode, previousNode, rt);
            else if(nReferrals == 3)
                nextReferralVec = refer_next_node<3>(adjlist, currentNode, previousNode, rt);
            else RASSERT(0)

            for(int i=0; i < nextReferralVec.size(); ++i){
                const auto nextReferral    = nextReferralVec[i];
                const int  nextNode        = get<0>(nextReferral);
                toBeVisited.push( QObject( nextNode, currentNode, get<1>(nextReferral), wave+1) );
            }
        }
    }

    RASSERT(counter>0);
    if(counter != nSamples){
        Rcpp::Rcout<<"in counter less than nSamples"<<std::endl;
        Matrix newlogs = Matrix(counter, logs.ncol());
        for(int i=0; i<counter; ++i)
            for(int j=0; j<logs.ncol(); ++j)
                newlogs(i,j) = logs(i,j);
        logs = newlogs;
    }
    return logs;
}


// [[Rcpp::export]]
std::vector<std::vector<double> >  adj2list(SEXP X_) 
{

    typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
    const MapMatd A(Rcpp::as<MapMatd>(X_));

    RASSERT( A.cols() == A.rows() );
    AdjList adjlist;	
    for (int i=0; i<A.rows(); ++i){
        vector<double> vec;
        for (int j=0; j<A.cols(); ++j)
            if ( A(i,j) != 0 )
                vec.push_back( j );
        adjlist.push_back( vec );
    }

    for (auto list:adjlist)
        sort(list.begin(), list.end());
    return adjlist;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rdssim_cpp(Rcpp::List rcpp_adjlist, std::string rType, 
        bool Markovian, 
        int nSamples, int nReferrals, int seedNode, int rseed) {

    //NOTE: set the random seed
    srand(rseed);

    //NOTE: nReferrals == 1; chain
    //NOTE: nReferrals >  1; a tree
    
    AdjList adjlist;
    for(auto l:rcpp_adjlist)
        adjlist.push_back( Rcpp::as<AdjList::value_type>(l) );

    //FIXME add it as an input argument
    bool sortAdjacancyList = false;
    if( sortAdjacancyList )
        for (auto l:adjlist)
            sort(l.begin(), l.end());

    RASSERT( nSamples   > 0 );
    RASSERT( nReferrals > 0 );
    //NOTE: seedNode; It has to be modified from R array format that indexing starts from 1 to cpp format
    seedNode = seedNode - 1;
    if( seedNode < 0 || seedNode >= adjlist.size() )
        Rf_error("The seed node is not in the range!");

    ReferralType rt;
    if(rType == "sRW")
        rt = SIMP_RW;
    else if(rType == "acRW")
        rt = AC_RW;
    else 
        Rf_error("Undefined referral type!");

    Matrix logMatrix(nSamples, 4);
    colnames(logMatrix) = Rcpp::CharacterVector::create("participant_i", "participant_i+1", "weight", "wave");
    if(nReferrals == 1)
        sim_referral_chain(adjlist, nSamples, seedNode, rt, rseed, logMatrix);
    else
        sim_referral_tree(adjlist, nSamples, nReferrals, seedNode, rt, Markovian, rseed, logMatrix);

    return logMatrix;
}


