#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace fs = std::filesystem;


struct PopulationMetrics {
    double min{};
    double max{};
    double avg{};
    double std{};
    std::size_t size{}; // population size
};

class StatsInterface {
public:
    static const std::string statsPath;

    static PopulationMetrics calculate_population_metrics(const std::vector<double>& data) ;
    static bool create_directories_if_not_exists(const std::string& directoryPath);
    static void stats_for_multiple_trials(const std::string& filePath, const std::vector<double>& data); // open a file, save the statistical info, and then close it

    virtual ~StatsInterface() = default;
};
