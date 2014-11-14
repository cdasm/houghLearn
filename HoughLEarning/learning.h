#pragma once


#include <vector>
#include <unordered_map>
#include <algorithm>
#include <assert.h>
#include "basicsettings.inl"

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace HoughLearning;


using namespace std;

using namespace Eigen;


auto learningMapForInstance(const vector<vector<int> >& votes,const vector<int>& positions)->vector<unordered_map<int,vector< pair<pair<int,int>,int>>>>;