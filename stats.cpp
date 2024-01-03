#include "stats.h"

using namespace std;


const std::string StatsInterface::statsPath = "stats";

PopulationMetrics StatsInterface::calculate_population_metrics(const std::vector<double> &data) {
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

bool StatsInterface::create_directories_if_not_exists(const string &directoryPath) {
    if (!fs::exists(directoryPath)) {
        try {
            fs::create_directories(directoryPath);
//            std::cout << "Directory created successfully: " << directoryPath << std::endl;
            return true;
        } catch (const std::exception& e) {
//            std::cerr << "Error creating directory: " << e.what() << std::endl;
            return false;
        }
    } else {
//        std::cout << "Directory already exists: " << directoryPath << std::endl;
        return true;
    }
}

void StatsInterface::stats_for_multiple_trials(const std::string& filePath, const std::vector<double>& data) {
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
