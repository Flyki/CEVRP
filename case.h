#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <float.h>
#include <string>
#include <stdio.h>

using namespace std;

class Case {
public:
	Case(string, int);
	~Case();
	double getDistance(int, int);
	double getEnergyDemand(int, int);
	double calculateRouteDistance(vector<int>);
	double calculateRouteDistance(int*, int);
	double calculateRouteDistance(int*);
	int findNearestStation(int);
	int findNearestStation(int, int);
	int findNearestStationFeasible(int, int, double);
	int findNearestStationFeasible2(int, int, double);
	void writeAllPositions();
	void drawARoute(vector<int>, string);
	void testTheStationReach();
	void checkASoluton(string);
	vector<int> findTheNonDominatedStations(int, int);
	vector<set<int>> getCandiList(int candino);
	vector<vector<int>> getCandiList2(int candino);

	int depotNumber;
	int customerNumber;
	int stationNumber;
	int vehicleNumber;
	int depot;
	vector<pair<double, double>> positions;
	double** distances;
	vector<int>** bestStations;
	int** bestStation;
	vector<double> demand;
	double maxC;
	double maxQ;
	double conR;
	double maxDis;
	double totalDem;
	string filename;
	bool posflag;
	int ID;
	vector<vector<int>> candidatelist;
};
