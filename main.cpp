//
// Created by Yinghao Qin on 07/12/2023.
//

#include <iostream>
#include <fstream>
#include <thread>
#include <sstream>

#include "case.h"
#include "utilities.h"
#include "ant.h"
#include "CACO.h"
#include "struct.h"
#include "BACO2.h"

#define MAX_TRIALS 20

using namespace std;

const string DATA_PATH = "instances/";

struct PopulationMetrics {
    double min{};
    double max{};
    double avg{};
    double std{};
    std::size_t size{}; // population size
};


PopulationMetrics calculate_population_metrics(const std::vector<double> &data) {
    PopulationMetrics metrics;

    metrics.size = data.size();

    if (data.empty()) {
        metrics.min = 0.0;
        metrics.max = 0.0;
        metrics.avg = 0.0;
        metrics.std = 0.0;
    } else {
        // Calculate min, max, mean, and standard deviation for feasible data
        metrics.min = *std::min_element(data.begin(), data.end());
        metrics.max = *std::max_element(data.begin(), data.end());

        double sum = 0.0;
        for (double value : data) {
            sum += value;
        }
        metrics.avg = sum / static_cast<double>(data.size());

        double sumSquaredDiff = 0.0;
        for (double value : data) {
            double diff = value - metrics.avg;
            sumSquaredDiff += diff * diff;
        }
        metrics.std = (data.size() == 1) ? 0.0 : std::sqrt(sumSquaredDiff / static_cast<double>(data.size() - 1));
    }

    return metrics;
}

void stats_for_multiple_trials(const std::string& filePath, const std::vector<double>& data) {
    std::ofstream logStats;

    logStats.open(filePath);

    std::ostringstream oss;
    for (auto& perf:data) {
        oss << fixed << setprecision(2) << perf << endl;
    }
    PopulationMetrics metric = calculate_population_metrics(data);
    oss << "Mean " << metric.avg << "\t \tStd Dev " << metric.std << "\t " << endl;
    oss << "Min: " << metric.min << "\t " << endl;
    oss << "Max: " << metric.max << "\t " << endl;
    logStats << oss.str() << flush;

    logStats.close();
}



int main(int argc, char *argv[]) {
    int run;

    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <problem_instance_filename> <pop_init> <confidence_based_selection> <representation>" << endl;
        return 1;
    }

    string instanceName(argv[1]);
    int isCan = std::stoi(argv[2]);; // population initialization, 1 for "buildSolutionFromCandidates", otherwise "buildSolution"
    int isRA = std::stoi(argv[3]); // weather using confidence-based selection strategy, 1 means yes, 0 means no
    int representation = std::stoi(argv[4]); // 1 represents order-split, 2 represents direct with local search

    std::vector<double> perfOfTrials(MAX_TRIALS);
    std::vector<std::thread> threads;

    // Define a function to perform the threaded work
    auto thread_function = [&](int run) {
        Case* instance = new Case(DATA_PATH + instanceName, run);

        double timer; // the stop criteria of the algorithm max execution time,
        if (instance->customerNumber <= 100) {
            timer = 1.0 /100;
        } else if (instance->customerNumber <= 915) {
            timer = 2.0 /100;
        } else {
            timer = 3.0 /100;
        }
        double afr = 0.8; // confidence ratio of recharging gammaR, 0.8

        CACO* caco = new CACO(instance, run, isCan, isRA, representation, timer, afr);
        caco->run();

        perfOfTrials[run - 1] = caco->bestSolution->fit;

        delete caco;
        delete instance;
    };

    // Launch threads
    for (run = 1; run <= MAX_TRIALS; ++run) {
        threads.emplace_back(thread_function, run);
    }

    // Wait for threads to finish
    for (auto& thread : threads) {
        thread.join();
    }


    string instancePredix = instanceName.substr(0, instanceName.find_last_of('.'));
    string directoryPath = STATS_PATH + "/" + instancePredix;
    string filepath = directoryPath + "/" + "stats." + instancePredix + ".txt";
    stats_for_multiple_trials(filepath, perfOfTrials);

    return 0;
}