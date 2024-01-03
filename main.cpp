//
// Created by Yinghao Qin on 07/12/2023.
//

#include <iostream>
#include <thread>

#include "case.h"
#include "CACO.h"
#include "BACO2.h"
#include "stats.h"

#define MAX_TRIALS 20

using namespace std;

// Enum for Algorithms
enum class Algorithm {BACO, CBACO};

const string DATA_PATH = "instances/";
const string STATS_PATH = "stats";


int main(int argc, char* argv[]) {
    int run;
    if (argc != 3 && argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <algorithm: baco, cbaco> <problem_instance_filename>\n";
        std::cerr << "       " << argv[0] << " <algorithm: baco, cbaco> <problem_instance_filename> <pop_init> <confidence_based_selection> <representation>\n";
        return 1; // Exit with an error code
    }

    std::string algorithmStr(argv[1]);
    Algorithm algorithm;

    // Convert algorithm string to enum
    if (algorithmStr == "baco") {
        algorithm = Algorithm::BACO;
    } else if (algorithmStr == "cbaco") {
        algorithm = Algorithm::CBACO;
    } else {
        std::cerr << "Error: Unknown algorithm '" << algorithmStr << "'\n";
        return 1; // Exit with an error code
    }

    std::string instanceName(argv[2]);
    std::vector<double> perfOfTrials(MAX_TRIALS);
    std::vector<std::thread> threads;

    switch (algorithm) {
        case Algorithm::BACO: {
            auto thread_function_1 = [&](int run) {
                Case* instance = new Case(DATA_PATH + instanceName, run);
                auto* baco = new BACO2(instance, run);
                baco->run();

                perfOfTrials[run - 1] = baco->gbestf;

                delete instance;
                delete baco;
            };
            // Launch threads
            for (run = 2; run <= MAX_TRIALS; ++run) {
                threads.emplace_back(thread_function_1, run);
            }
            thread_function_1(1);

            // Wait for threads to finish
            for (auto& thread : threads) {
                thread.join();
            }
            break;
        }

        case Algorithm::CBACO: {
            int isCan = std::stoi(argv[3]); // population initialization, 1 for "buildSolutionFromCandidates", otherwise "buildSolution"
            int isRA = std::stoi(argv[4]); // weather using confidence-based selection strategy, 1 means yes, 0 means no
            int representation = std::stoi(argv[5]); // 1 represents order-split, 2 represents direct with local search

            auto thread_function_2 = [&](int run) {
                Case* instance = new Case(DATA_PATH + instanceName, run);

                int timer;
                if (instance->customerNumber <= 100) {
                    timer = 1 * 36;
                } else if (instance->customerNumber <= 915) {
                    timer = 2 * 36;
                } else {
                    timer = 3 * 36;
                }
                double afr = 0.8; // confidence ratio of recharging gammaR, 0.8

                CACO* caco = new CACO(instance, run, isCan, isRA, representation, timer, afr);
                caco->run();

                perfOfTrials[run - 1] = caco->bestSolution->fit;

                delete caco;
                delete instance;
            };

            // Launch threads
            for (run = 2; run <= MAX_TRIALS; ++run) {
                threads.emplace_back(thread_function_2, run);
            }
            thread_function_2(1);

            // Wait for threads to finish
            for (auto& thread : threads) {
                thread.join();
            }
            break;
        }
    }

    string instancePrefix = instanceName.substr(0, instanceName.find_last_of('.'));
    string directoryPath;
    if (algorithm == Algorithm::CBACO) {
        int representation = std::stoi(argv[5]);
        if (representation == 1) {
            directoryPath = STATS_PATH + "/" + algorithmStr + "-i" + "/" + instancePrefix; // Use algorithmStr
        } else {
            directoryPath = STATS_PATH + "/" + algorithmStr + "-d" + "/" + instancePrefix; // Use algorithmStr
        }
    } else {
        directoryPath = STATS_PATH + "/" + algorithmStr + "/" + instancePrefix;
    }
    string filepath = directoryPath + "/" + "stats." + instancePrefix + ".txt";
    StatsInterface::stats_for_multiple_trials(filepath, perfOfTrials);

    return 0;
}