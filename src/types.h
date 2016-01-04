#ifndef TYPES_H
#define TYPES_H

#define RASSERT(condition){if(!(condition)){Rf_error(("error@"+to_string(__LINE__)+":"+__FILE__).c_str());}}

#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include <cstring>
#include <random>
#include <ctime>
#include <queue>
#include <set>

//Rcpp headers
#include <Rcpp.h>
#include <RcppEigen.h>

using std::cout; 
using std::endl;
using std::min; 
using std::vector;
using std::string;
using std::to_string;
using std::range_error;
using std::for_each;
using std::sort;
using std::pair;
using std::tuple;
using std::make_tuple;
using std::get;
using std::queue;

using Rcpp::Rcout;

typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Rcpp::NumericMatrix Matrix;

typedef tuple<int,double> Referral;
typedef vector<double> LogVector;
typedef vector<vector<double> > LogVectors;
typedef vector<int> Chain;
typedef vector<string> NameVector;
typedef vector<vector<double> > AdjList;
typedef vector<vector<double> > AdjWeightList;
typedef vector<double> TraitList;
typedef std::pair<int,int> NodePair;
typedef std::vector<NodePair> NodePairVector;

typedef std::default_random_engine RandEngine;
typedef std::uniform_int_distribution<int> IntDist;

#endif
