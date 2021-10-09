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
using namespace std;

#define INFEASIBLE 1000000000

vector<int> convexHull(vector<int> nodesemp, Case* instance);

//insert every time really needs to charge
pair<vector<int>, double> insertStationByNecessaryCharge(vector<int> route, Case* instance);

//insert a station on each link and remove the one cause max extra distance
pair<vector<int>, double> insertStationByRemove(vector<int> route, Case* instance);

//regular 2opt
bool opt2noStation(vector<int>& route, Case* instance);
bool opt2noStation2(int* route, int length, Case* instance);
bool orNoStation(int* route, int length, Case* instance);
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
