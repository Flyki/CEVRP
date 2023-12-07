//
// Created by Yinghao Qin on 07/12/2023.
//

#include <iostream>
#include <fstream>

#include "case.h"
#include "utilities.h"
#include "ant.h"
#include "CACO.h"
#include "struct.h"


using namespace std;

const string DATA_PATH = "instances/";

int main(int argc, char *argv[]) {

    string instanceName(argv[1]);

    int ID = 1;
    Case* instance = new Case(DATA_PATH + instanceName, ID);

    int seed = 1;
    int isCan = 1; // population initialization, 1 for "buildSolutionFromCandidates", otherwise "buildSolution"
    int isRA = 1; // weather using confidence-based selection strategy, 1 means yes, 0 means no
    int representation = 1; // 1 represents order-split, 2 represents direct with local search
    double timer; // the stop criteria of the algorithm max execution time,
    if (instance->customerNumber <= 100) {
        timer = 1.0 /100;
    } else if (instance->customerNumber <= 915) {
        timer = 2.0 /100;
    } else {
        timer = 3.0 /100;
    }

    double afr = 0.8; // confidence ratio of recharging gammaR, 0.8

    CACO* caco = new CACO(instance, seed, isCan, isRA, representation, timer, afr);
    caco->run();


    delete instance;
    delete caco;
    return 0;
}

