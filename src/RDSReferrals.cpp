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


enum ReferralType  {SIMP_RW, AC_RW};
enum AC_Type       {ONE, TWO, BOTH};

// returns a random integer X in the set {l, l+1, ..., u}.
int random_int(const int l, const int u){
	RASSERT( l <= u );
	double uR = l + (u+1 - l)*(double) rand()/((double) RAND_MAX + 1);
	return (int) uR;
}


inline bool are_not_neighbor(const AdjList &adjList, const int node1, const int node2){

	if(node1 == node2) 
		return 0;

	const auto l = adjList[node1];
	return !(std::binary_search(l.begin(), l.end(), node2));
}

inline bool are_neighbor(const AdjList &adjList, const int node1, const int node2){

	if(node1 == node2) 
		return 1;

	const auto l = adjList[node1];
	return std::binary_search(l.begin(), l.end(), node2);
}

void refer_next_node(const AdjList &, const AdjList &, const int, const int, const int, const ReferralType, vector<Referral>&);


void enumerate_ac_rw_possible_referrals(const AdjList& adjList, const int cNode, 
		vector<int>& possibleReferrals, const AC_Type acType){

	const auto list = adjList[cNode];
	if( acType == BOTH || acType == ONE ){
		for(auto item1:list)
			for(auto item2:list)
				if ( are_not_neighbor(adjList, item1, item2) )
					possibleReferrals.push_back( item1 );
	}

	if( acType == BOTH || acType == TWO ){
		for(auto item1:list)
			for(auto item2:adjList[item1])
				if ( are_not_neighbor(adjList, item2, cNode) )
					possibleReferrals.push_back( item1 );
	}
}

void refer_next_node_ac_rw(const AdjList &adjList, const AdjList &acAdjList, 
		const int cNode, const int prevNode, 
		const int nReferrals, const ReferralType rt, 
		vector<Referral>& refVec
		){


	//vector<int> possibleReferrals;
	//possibleReferrals.reserve(1024);
	//enumerate_ac_rw_possible_referrals(adjList, cNode, possibleReferrals);

	const auto& possibleReferrals = acAdjList[cNode];
	const int nCandids = possibleReferrals.size() ;
	const double weight = 1;

	if( !nCandids ){//Handling special cases
		Rf_warning(("cN "+to_string(cNode)+". Zero candid list!").c_str());
		refer_next_node(adjList, acAdjList, cNode, prevNode, nReferrals, (ReferralType) 0, refVec);
		return;
	}

	for(int i=0; i < nReferrals; i++){
		const int id = random_int(0, nCandids-1);
		RASSERT(id < possibleReferrals.size() );
		const int nextNode = possibleReferrals[id];
		refVec.push_back( make_tuple(nextNode, weight) );
	}
}


void refer_next_node(const AdjList &adjList, 
		const AdjList &acAdjList,
		const int cNode, const int prevNode, 
		const int nReferrals, const ReferralType rt,
		 vector<Referral>& refVec){

	RASSERT( cNode < adjList.size() && cNode >= 0);
	RASSERT( prevNode < adjList.size() && prevNode >= 0);
	RASSERT( nReferrals > 0 );
	//RASSERT( nReferrals < 4 );

	if( rt == AC_RW ) {
		refer_next_node_ac_rw (adjList, acAdjList, cNode, prevNode, nReferrals, rt, refVec);
		return;
	}


	const auto list      = adjList[cNode];
	const int nNeighbors = list.size();
	const double weight  = 1/(double) nNeighbors;

	for(int i=0; i<nReferrals; i++){
		const int firstIdx = random_int(0, nNeighbors-1);
		const int firstNode = list[ firstIdx ];
		refVec.push_back( make_tuple(firstNode, weight) );
	}

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

Matrix sim_referral_chain(const AdjList &adjList,
		const AdjList &acAdjList,
		const ReferralType rt, 
		const int chainLength, 
		const int seedNode, const unsigned int rseed, 
		Matrix &logs){


	RASSERT(logs.nrow() == chainLength);
	RASSERT(logs.ncol() == 4);

	if(adjList[seedNode].size() == 0)
		Rf_error("The seed node is an isolated node!");

	int currentNode = seedNode, previousNode = seedNode;

	for(int level = 0; level<chainLength; ++level){

		const int nNeighbors = adjList[currentNode].size();
		if (nNeighbors == 0)
			Rf_error("Zero neighbors! That should not happen.");

		vector<Referral> nextReferralVec;
	       	refer_next_node(adjList, acAdjList, currentNode, previousNode, 1, rt, nextReferralVec);
		const auto nextReferral  = nextReferralVec[0];

		//NOTE: I assume that the node that enters the sampling is referred directly by a participant.
		//RASSERT( are_neighbor(adjList, currentNode, get<0>( nextReferral ) ) );

		fill_log(logs, currentNode, nextReferral, level, level); 

		previousNode = currentNode;
		currentNode  = get<0>(nextReferral);
	}
	return logs;
}


Matrix sim_referral_tree( const AdjList &adjList,
		const AdjList &acAdjList,
		const ReferralType rt, const bool wReplacement,
		const int nSamples,    int nReferrals,
		const int seedNode,    const unsigned int rseed, 
		Matrix &logs
		){

	if (adjList[seedNode].size() == 0)
		Rf_error("The seed node is an isolated node!");

	if (acAdjList[seedNode].size() == 0)
		Rf_error("The seed node is an isolated node in AcAdjList!");

	RASSERT(logs.nrow() == nSamples);
	RASSERT(logs.ncol() == 4);

	bool Markovian = true;
	if(nReferrals < 0){
		nReferrals = -1*nReferrals;
		Markovian  = false;
	}

	int currentNode = seedNode, previousNode = seedNode, wave=0;
	double weight = 1;

	std::set<int> visitedNodes;
	visitedNodes.insert(currentNode);

	//node, parent, weight, wave number
	typedef tuple<int, int, double, int> QObject;
	queue<QObject> toBeVisited;
	toBeVisited.push( QObject(currentNode, previousNode, weight, wave) );

	int counter=0;
	for(; counter<nSamples; ++counter){

		RASSERT( adjList[currentNode].size()   > 0 );
		RASSERT( acAdjList[currentNode].size() > 0 );

		if( toBeVisited.size() == 0 )
			break;

		const auto qObj = toBeVisited.front(); 
		toBeVisited.pop();

		currentNode  = get<0>(qObj);
		previousNode = get<1>(qObj);
		weight       = get<2>(qObj);
		wave         = get<3>(qObj);

		//RASSERT( ! are_not_neighbor(adjList, currentNode, previousNode) );


		//if(!wReplacement){
		//    if( visitedNodes.find(currentNode) != visitedNodes.end() ){
		//        --counter;
		//        continue;
		//    }
		//    visitedNodes.insert(currentNode);
		//}

		fill_log(logs, previousNode, make_tuple(currentNode, weight), wave, counter); 

		vector<Referral> nextReferralVec; 
		if( wReplacement ){
			refer_next_node(adjList, acAdjList, currentNode, previousNode, nReferrals, rt, nextReferralVec);
			for(int i=0; i < nReferrals; ++i){
				const auto nextReferral    = nextReferralVec[i];
				const auto nextNode        = get<0>(nextReferral);
				toBeVisited.push( QObject( nextNode, currentNode, get<1>( nextReferral ), wave+1) );
			}
		} else {
			const int maxNReferrals = nReferrals; //10*nReferrals;
			refer_next_node(adjList, acAdjList, currentNode, previousNode, maxNReferrals, rt, nextReferralVec);
			for(int i=0, refIdx=0; i < nReferrals && refIdx < maxNReferrals; ++i, ++refIdx){

				auto nextReferral = nextReferralVec[refIdx];
				auto nextNode     = get<0>(nextReferral);

				//while( visitedNodes.find(nextNode) != visitedNodes.end() 
				//		&& refIdx < (maxNReferrals -1) ){
				//	nextReferral = nextReferralVec[++refIdx];
				//	nextNode     = get<0>(nextReferral);
				//}

				if( visitedNodes.find(nextNode) == visitedNodes.end() ){
					toBeVisited.push( QObject( nextNode, currentNode, get<1>( nextReferral ), wave+1) );
					visitedNodes.insert(nextNode);
				}
			}
		}

	}

	RASSERT(counter>0);
	if(counter != nSamples){
		//Rcpp::Rcout<<"in counter less than nSamples"<<std::endl;
		Matrix newlogs = Matrix(counter, logs.ncol());
		for(int i=0; i<counter; ++i)
			for(int j=0; j<logs.ncol(); ++j)
				newlogs(i,j) = logs(i,j);
		logs = newlogs;
	}
	return logs;
}

// [[Rcpp::export]]
Rcpp::List adj2list(SEXP X_, std::string acType){

	typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
	const MapMatd A(Rcpp::as<MapMatd>(X_));

	RASSERT( A.cols() == A.rows() );
	AdjList adjList;	
	for (int i=0; i<A.rows(); ++i){
		vector<int> vec;
		for (int j=0; j<A.cols(); ++j)
			if ( A(i,j) != 0 )
				vec.push_back( j );
		adjList.push_back( vec );
	}

	for (auto list:adjList)
		sort(list.begin(), list.end());

	AC_Type act;
	if( acType == "One")
		act = ONE;
	else if( acType == "Two")
		act = TWO;
	else if( acType == "Both")
	      act = BOTH;	
	else
		Rf_error("Undefined anti clustering referral type!");
	
	AdjList acAdjList;	
	for (int node=0; node<adjList.size(); node++){
		vector<int> vec;
		enumerate_ac_rw_possible_referrals(adjList, node, vec, act);
		acAdjList.push_back( vec );
	}


	Rcpp::List L;
	L = Rcpp::List::create( Rcpp::Named("AdjList") = adjList,  
				Rcpp::Named("AcAdjList") = acAdjList  );
	return L;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rdssim_cpp(Rcpp::List rcpp_adjList, 
		Rcpp::List rcpp_acAdjList, 
		std::string referralType, 
		bool wReplacement, 
		int nSamples, int nReferrals, int seedNode, int rseed) {

	//NOTE: set the random seed
	srand(rseed);


	ReferralType rt;
	if(referralType == "sRW")
		rt = SIMP_RW;
	else if(referralType == "acRW")
		rt = AC_RW;
	else 
		Rf_error("Undefined referral type!");

	//NOTE: nReferrals == 1; chain
	//NOTE: nReferrals >  1; a tree

	AdjList adjList, acAdjList;
	for(auto l:rcpp_adjList)
		//adjList.push_back( Rcpp::as<AdjList::value_type>(l) );
		adjList.push_back( Rcpp::as<vector<int> >(l) );

	for(auto l:rcpp_acAdjList)
		acAdjList.push_back( Rcpp::as<AdjList::value_type>(l) );


	//FIXME add it as an input argument
	bool sortAdjacancyList = false;
	if( sortAdjacancyList )
		for (auto l:adjList)
			sort(l.begin(), l.end());

	RASSERT( nSamples   > 0 );
	RASSERT( nReferrals > -4 );

	if(nReferrals == 1 && wReplacement==false)
		nReferrals = -1; // it is implemented in sim_referral_tree 

	//NOTE: seedNode; It has to be modified from R array format that indexing starts from 1 to cpp format
	seedNode = seedNode - 1;
	if( seedNode < 0 || seedNode >= adjList.size() )
		Rf_error("The seed node is not in the range!");


	Matrix logMatrix(nSamples, 4);
	colnames(logMatrix) = Rcpp::CharacterVector::create("participant_i", "participant_i+1", "weight", "wave");
	if(nReferrals == 1)
		sim_referral_chain (adjList, acAdjList, rt, nSamples, seedNode, rseed, logMatrix);
	else
		sim_referral_tree  (adjList, acAdjList, rt, wReplacement, nSamples, nReferrals, seedNode, rseed, logMatrix);

	return logMatrix;
}


