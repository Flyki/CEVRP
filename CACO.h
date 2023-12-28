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
#include <filesystem>

#include "case.h"
#include "utilities.h"
#include "ant.h"

const string STATS_PATH = "stats/baco";

namespace fs = std::filesystem;
using namespace std;

class CACO {
public:
    static bool create_directories_if_not_exists(const std::string& directoryPath); // by Yinghao

    //for representation 1 represents order-split, 2 represents direct with local search
    CACO(Case* instance, int seed, int isCan, int isRA, int representation, double timer, double afr);
    ~CACO();
    void run();
    //traditional representations
    void buildSolutions();
    //order split representations
    void buildSolutions2();
    ////traditional representations
    void buildSolutionsFromCandi();
    //order split representations
    void buildSolutionsFromCandi2();
    //traditional representations
    void evaluateAll();
    //for ordersplit
    void evaluateAll2();
    //traditional representations
    void evaluateSome();
    //for ordersplit only use 2opt
    void evaluateSome2();

    void evaluateSomeForOnlyLocalSearch1();
    void evaluateSomeForOnlyLocalSearch0();
    void evaluateSomeForOnlyFixing1();
    void evaluateSomeForOnlyFixing0();
    
    void generateASolutionGreedy(Ant*);
    double fixOneSolution(Ant*);

    Case* instance;
    ofstream result;
    ofstream sofile;
    int antno;
    vector<Ant*> ants;
    Ant* bestSolution;
    double** pher;
    int cdnumber;
    double maxt;
    double mint;
    double t0;
    bool isCan;
    bool isRA;
    bool representation;
    int ibest;
    double afr;

    default_random_engine gen;
    uniform_real_distribution<double> udis;
    normal_distribution<double> ndis;
    clock_t staTime;
    clock_t endTime;
    vector<vector<int>> candidatelist;
    vector<double> refiningImprovement;
    double maxRefine;
    double minRepair;
    long int usedFes;
    long int refined;
    long int repaired;
    double timerate;
};
