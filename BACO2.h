//always start from zero.
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <random>
//#include <sys/time.h>
#include <time.h>
#include <queue>
#include "case.h"
#include "utilities.h"

using namespace std;

class BACO2 {
public:
	BACO2(Case*, int);
	~BACO2();
	void run();
	void buildSolutionsByCL();
	void buildSolutionsByAll();
	void evaluateAndUpdatePher();
	pair<vector<int>, double> interpretACircleVec(int*);

	Case* instance;
	ofstream result;
	ofstream sofile;
	int ngen;
	int antno;
	int** ants;
	double* fit;
	int* gbest;
	double gbestf;
	double t0;
	double** pher;
	int cdnumber;
	double maxt;
	double mint;
	int ibest;
	vector<int> bestSolution;
	default_random_engine gen;
	uniform_real_distribution<double> udis;
	long usedFes;
	int candinumber;
	//timeval t1, t2;
	clock_t t1, t2;
	vector<vector<int>> candidatelist;
};