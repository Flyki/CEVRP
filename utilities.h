#pragma once
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_set>
#include <list>
#include <sstream>
#include "case.h"
#include "struct.h"
#include "ant.h"

using namespace std;

#define INFEASIBLE 1000000000
#define BIGVALUE 1000000

double getLargest(vector<double> x);
double getSmallest(vector<double> x);
vector<int> convexHull(vector<int> nodesemp, Case* instance);

//insert every time really needs to charge
pair<vector<int>, double> insertStationByNecessaryCharge(vector<int> route, Case* instance);

//insert a station on each link and remove the one cause max extra distance
pair<vector<int>, double> insertStationByRemove(vector<int> route, Case* instance);
pair<vector<int>, double> insertStationByRemove2(vector<int> route, Case* instance);
//regular 2opt
bool opt2noStation(vector<int>& route, Case* instance);
bool opt2noStation2(int* route, int length, double& fitv, Case* instance);
bool orNoStation(int* route, int length, double& fitv, Case* instance);
bool exchangeNoStation(vector<int>& route, Case* instance);
//include start and end
void reverseAnArray(int*, int start, int end);
void moveItoJ(int*, int i, int j);

pair<vector<int>, double> greedymethod(Case* instance);

//true if equal, false if not equal
bool compareArr(vector<int>, int*);

//prins split
vector<vector<int>> prinsSplit(vector<int> x, Case* instance);

//calcualte dificit and distance, prepare for insertStationByMinDeficit
//route should start and end 0
pair<double, double> deficitAndDistance(vector<int>& route, Case* instance);
//choose the one reduce max deficit and cause min extra distance
pair<vector<int>, double> insertStationByMinDeficit(vector<int> route, Case* instance);

pair<vector<int>, double> insertStationByRemove3(vector<int> route, Case* instance);

//pair<vector<int>, double> insertStationByEnumerate(vector<int> route, Case* instance);
pair<vector<int>, double> insertStationByEnumeration(vector<int> route, Case* instance);
pair<vector<int>, double> insertStationByEnumeration(vector<int> route, Case* instance, int ub);
void tryACertainN(int mlen, int nlen, int* chosenSta, int* chosenPos, vector<int>& finalRoute, double& finalfit, int curub, vector<int>& route, vector<double>& accumulateDis, Case* instance);
//void evaluateAnInsertion(int* chosenPos, int curub, vector<int>& finalRoute, double& finalfit, vector<int>& route, Case* instance);
pair<vector<int>, double> insertStationByEnumerationWithStation(vector<int> route, Case* instance);
pair<vector<int>, double> insertStationBySimpleEnumeration(vector<int> route, Case* instance);
void tryACertainN(int mlen, int nlen, int* chosenPos, vector<int>& finalRoute, double& finalfit, int curub, vector<int>& route, vector<double>& accumulateDis, Case* instance);
double insertStationBySimpleEnumeration2(vector<int> route, Case* instance);
//starts and ends at depot
bool opt2starNoStation(Ant* ant, Case* instance);
//maybe faster
struct pair_hash
{
    size_t operator() (pair<int, int> const & apair) const {
        return apair.first * 256 + apair.second;
    }
};

bool opt2starNoStation2(Ant* ant, Case* instance);
bool opt2starNeighborNoStation(Solution route, Case* instance);
void opt2ToAnt(Ant* ant, Case* instance);
void orToAnt(Ant*, Case*);

double insertStationBySimpleEnumerationArray(int* route, int length, Case* instance);
void tryACertainNArray(int mlen, int nlen, int* chosenPos, double& finalfit, int curub, int* route, int length, vector<double>& accumulatedDis, Case* instance);
double insertStationByRemoveArray(int* route, int length, Case* instance);
void prinsSplitAnt(Ant* x, Case* instance);

//include 2opt with station, 2opt* with station, or with station
void finalRefine(Ant* x, Case* instance);
bool opt2withStation(int* route, int length, Case* instance);
