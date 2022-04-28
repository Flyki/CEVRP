#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <random>
#include <time.h>
#include "case.h"
#include "utilities.h"

using namespace std;

class Ant {
public:
	Ant(int routeCap, int nodeCap) {
		this->routeCap = routeCap;
		this->nodeCap = nodeCap;
		this->route = new int* [routeCap];
        this->demsum = new double [routeCap];
        memset(demsum, 0, sizeof(double) * routeCap);
		this->circle = new int[nodeCap];
		for (int i = 0; i < routeCap; i++) {
			this->route[i] = new int[nodeCap];
			memset(this->route[i], 0, sizeof(int) * nodeCap);
		}
		this->nodeNum = new int[routeCap];
		memset(this->nodeNum, 0, sizeof(int) * routeCap);
		routeNum = 0;
		this->fit = 0;
	}
	~Ant() {
		for (int i = 0; i < routeCap; i++) {
			delete[] this->route[i];
		}
        delete[] this->demsum;
		delete[] this->route;
		delete[] this->nodeNum;
		delete[] this->circle;
	}
	void reset() {
		memset(this->nodeNum, 0, sizeof(int) * routeCap);
		memset(this->demsum, 0, sizeof(double) * routeCap);
		this->fit = 0;
		routeNum = 0;
	}
	void copyASolutionIntoAnt(vector<vector<int>>& as, double afit, double* dem) {
		this->fit = afit;
		this->routeNum = (int)as.size();
		for (int i = 0; i < (int)as.size(); i++) {
			this->nodeNum[i] = (int)as[i].size();
			for (int j = 0; j < (int)as[i].size(); j++) {
				this->route[i][j] = as[i][j];
			}
		}
		memcpy(this->demsum, dem, sizeof(double) * routeCap);
	}
	vector<vector<int>> getTheRoutes() {
		vector<vector<int>> allroutes(routeNum);
		for (int i = 0; i < routeNum; i++) {
			allroutes[i].resize(nodeNum[i]);
			for (int j = 0; j < nodeNum[i]; j++) {
				allroutes[i][j] = route[i][j];
			}
		}
		return allroutes;
	}
	/*Ant operator=(const Ant& anant) {
		this->fit = anant.fit;
		this->routeNum = anant.routeNum;
		memcpy(this->nodeNum, anant.nodeNum, sizeof(int) * this->routeNum);
		for (int i = 0; i < this->routeNum; i++) {
			memcpy(this->route[i], anant.route[i], sizeof(int) * this->nodeNum[i]);
		}
		return anant;
	}*/

	int** route;
	int* circle;
    double* demsum;
	int routeNum;
	int* nodeNum;
	int routeCap;
	int nodeCap;
	double fit;
};
