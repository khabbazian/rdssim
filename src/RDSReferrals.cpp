// Under GNU (>=2) license 
//
// flags for Rcpp 
// [[Rcpp::depends("RcppEigen")]]
// [[Rcpp::plugins(cpp11)]]

#include "types.h"

//NOTE: SIMP_RW:  A simple Markovian random walk.
//NOTE: AC_RW:  AC random walk choses uniformly from the set of nodes form wedge
//NOTE: with current node.

enum ReferralType{SIMP_RW, AC_RW};

// checks if nodes 1 is in the list of node2
bool is_wedge(const AdjList &adjlist, const int node1, const int node2){

    RASSERT( node1 < adjlist.size() );
    RASSERT( node2 < adjlist.size() );

    if(node1 == node2 ) 
        return 0;

    for(auto item:adjlist[node1])
        if(item == node2)
            return 0;
    return 1;
}


bool are_neighbor(const AdjList &adjlist, const int node1, const int node2){

    RASSERT( node1 < adjlist.size() );
    RASSERT( node2 < adjlist.size() );

    if(node1 == node2 ) 
        return 1;

    for(auto item:adjlist[node1])
        if( item == node2)
            return 1;
    return 0;
}

template<int nReferrals>
vector<Referral> refer_next_node(const AdjList &, const int, const int, const ReferralType, RandEngine &);
template<int nReferrals>
vector<Referral> refer_next_node_ac_rw(const AdjList &, const int, const int, const ReferralType, RandEngine &);


template<int nReferrals>
vector<Referral> refer_next_node_ac_rw(const AdjList &adjlist, const int cNode, const int prevNode, 
        const ReferralType rt,
        RandEngine &generator){

    int nextNode;
    double weight=1;
    const auto list = adjlist[cNode];

    NodePairVector possibleReferrals;
    int typeOneNCandids = 0, typeTwoNCandids = 0;
    for(auto item1:list){
        for(auto item2:list)
            if ( is_wedge(adjlist, item1, item2) )
                possibleReferrals.push_back( NodePair(item1, item2) );

        typeOneNCandids = possibleReferrals.size();

        for(auto item2:list)
            if ( is_wedge(adjlist, prevNode, item2) )
                possibleReferrals.push_back( NodePair(item2, prevNode) );

        typeTwoNCandids = possibleReferrals.size() - typeOneNCandids;
    }

    RASSERT( typeOneNCandids + typeTwoNCandids == possibleReferrals.size() );
    const int nCandids = typeOneNCandids + typeTwoNCandids;
    {
        int counter = 0;
        for(auto item:possibleReferrals)
            counter += get<0>(item) == nextNode ? 1 : 0;
        weight = counter/(double) nCandids;
    }

    //Handling special cases
    if( !nCandids ){
        Rf_warning(("cN "+to_string(cNode)+". Zero candid list!").c_str());
        return refer_next_node<nReferrals>(adjlist, cNode, prevNode, (ReferralType) 0, generator);
    }

    if( nReferrals == 1 || nCandids == 1 ){
        IntDist dist(0, nCandids-1);
        const int id = dist(generator);
        nextNode = get<0>( possibleReferrals[id] );
        return vector<Referral> { make_tuple(nextNode,weight) };
    } else if( nReferrals == 2 ){
        //@this moment I have no idea about the best scenario for this type of referral
        // here's just the first thing that came to mind
        IntDist distTypeOne(0, typeOneNCandids-1);
        IntDist distTypeTwo(0, typeTwoNCandids-1);
        const int idOne = distTypeOne(generator);
        const int idTwo = distTypeTwo(generator);

        const auto node1 = get<0>( possibleReferrals[idOne] );
        const auto node2 = get<0>( possibleReferrals[idTwo] );

        return vector<Referral> {make_tuple(node1,weight), make_tuple(node2,weight)};
    } else if( nReferrals == 3 ){
        vector<Referral> refVec;

        if( typeOneNCandids > 0 ){
            IntDist distTypeOne(0, typeOneNCandids-1);
            const int idOne = distTypeOne(generator);
            const auto node1 = get<0>( possibleReferrals[idOne] );
            const auto node2 = get<1>( possibleReferrals[idOne] );

            refVec.push_back( make_tuple(node1, weight) );
            refVec.push_back( make_tuple(node2, weight) );
        }

        if( typeTwoNCandids > 0 ){
            IntDist distTypeTwo(0, typeTwoNCandids-1);
            const int idTwo  = distTypeTwo(generator);
            const auto node3 = get<0>( possibleReferrals[idTwo] );
            refVec.push_back( make_tuple(node3, weight) );
        }

        return refVec;
    }
}


template<int nReferrals>
vector<Referral> refer_next_node(const AdjList &adjlist, const int cNode, const int prevNode, 
        const ReferralType rt,
        RandEngine &generator){

    RASSERT( cNode < adjlist.size() );
    RASSERT( prevNode < adjlist.size() );
    RASSERT( nReferrals < 4 && nReferrals > 0 );

    int nextNode;
    double weight=1;
    const auto list = adjlist[cNode];

    if( rt == SIMP_RW ) {
        const int nNeighbors = list.size();
        weight = 1/(double) nNeighbors;
        vector<Referral> refVec;
       
        //NOTE: scenario for referring the first node
        IntDist firstDist(0, nNeighbors-1);
        const int firstIdx = firstDist(generator); 
        nextNode  = list[ firstIdx ];
        RASSERT(nextNode < adjlist.size());
        refVec.push_back( make_tuple(nextNode,weight) );

        if (nReferrals == 1 || nNeighbors < 2) //NOTE: to generate Markov chains.
            return refVec;

        //NOTE: I avoided using for loop since I only consider nReferrals=2 and 3. 
        //Therefore it would be more readable in this way.

        //NOTE: scenario for referring the second node; a w/o replacement sampling from neighbors.
        IntDist secondDist(0, nNeighbors-2);
        int secondIdx = secondDist(generator);
        secondIdx = secondIdx < firstIdx ? secondIdx : secondIdx + 1;
        refVec.push_back( make_tuple( list[ secondIdx ], weight) );

        if (nReferrals == 2 || nNeighbors < 3) 
            return refVec;

        //NOTE: scenario for referring the third node; similar to the second one w/o replacement sampling form neighbors.
        IntDist thirdDist(0, nNeighbors-3);
        int thirdIdx = thirdDist(generator);

        thirdIdx  = thirdIdx < firstIdx  ? thirdIdx  : thirdIdx + 1;
        thirdIdx  = thirdIdx < secondIdx ? thirdIdx  : thirdIdx + 1;

        refVec.push_back( make_tuple(list[ thirdIdx ], weight) );
        return refVec;

    } else if( rt == AC_RW ) {
        return refer_next_node_ac_rw<nReferrals>(adjlist, cNode, prevNode, rt, generator);
    } else 
        RASSERT(0);
} 

inline void fill_log(Matrix &logs, const int currentNode, const Referral nextReferral, const int wave, const int l){
        //NOTE: currentNode; It has to be modified from R array format that indexing starts from 1
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

    if(adjlist[seedNode].size() == 0)
        Rf_error("The seed node is an isolated node!");

    RASSERT(logs.nrow() == chainLength);
    RASSERT(logs.ncol() == 4);

    std::default_random_engine generator (rseed);
    int currentNode = seedNode, previousNode = seedNode;

    for(int level = 0; level<chainLength; ++level){
        const int nNeighbors = adjlist[currentNode].size();

        if (nNeighbors == 0)
            Rf_error("Zero neighbors! That should not happen.");

        const auto nextReferralVec = refer_next_node<1>(adjlist, currentNode, previousNode, rt, generator);
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
        const unsigned int rseed, Matrix &logs){

    if (adjlist[seedNode].size() == 0)
        Rf_error("The seed node is an isolated node!");

    RASSERT(logs.nrow() == nSamples);
    RASSERT(logs.ncol() == 4);

    std::default_random_engine generator (rseed);
    int currentNode = seedNode, previousNode = seedNode, wave=0;
    double weight = 1;

    //node, parent, weight, wave number
    typedef tuple<int, int, double, int> QObject;
    queue<QObject> toBeVisited;
    toBeVisited.push( QObject(currentNode, previousNode, weight, wave) );

    int counter=0;
    for(int counter=0; counter<nSamples; ++counter){
        const int nNeighbors = adjlist[currentNode].size();
        RASSERT( nNeighbors > 0 );
        RASSERT( toBeVisited.size() > 0 );

        const auto qObj = toBeVisited.front(); 
        toBeVisited.pop();

        currentNode  = get<0>(qObj);
        previousNode = get<1>(qObj);
        weight       = get<2>(qObj);
        wave         = get<3>(qObj);

        fill_log(logs, previousNode, make_tuple(currentNode, weight), wave, counter); 

        for(int i=0; i < min(nReferrals, nNeighbors); ++i){
            const auto nextReferralVec = refer_next_node<1>(adjlist, currentNode, previousNode, rt, generator);
            const auto nextReferral    = nextReferralVec[0];
            toBeVisited.push( QObject( get<0>(nextReferral), currentNode, get<1>( nextReferral ), wave+1) );
        }
    }
    return logs;
}


// [[Rcpp::export]]
AdjList adj2list(SEXP X_) {

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

//Matrix rdssimMarkov(SEXP A_, string rType,
//        int nSamples, int nReferrals=1, int seedNode=1, int rseed=1)

// [[Rcpp::export]]
Matrix rdssimMarkov(Rcpp::List rcpp_adjlist, string rType,
        int nSamples, int nReferrals=1, int seedNode=1, int rseed=1)
{

    //NOTE: nReferrals == 1; chain
    //NOTE: nReferrals >  1; a tree

    RASSERT( nSamples   > 0 );
    RASSERT( nReferrals > 0 );
    RASSERT( seedNode   > 0 );

    //NOTE: seedNode; It has to be modified from R array format that indexing starts from 1
    seedNode = seedNode - 1;

    //AdjList adjlist = adj2list(A_);
    AdjList adjlist;
    for(auto l:rcpp_adjlist){
        AdjList::value_type tmp;
        for(auto item:Rcpp::NumericVector(l))
            tmp.push_back(item);
        adjlist.push_back( tmp );
    }

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
        sim_referral_tree(adjlist, nSamples, nReferrals, seedNode, rt, rseed, logMatrix);

    return logMatrix;
}


