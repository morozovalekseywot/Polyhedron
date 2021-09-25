//
// Created by 159-mrv on 6/21/18.
//

#include <utils/partitioning/splitter.h>
#include <utils/math.h>

void print_workload_info(const string &tabs, const vector<double>& workloads) {
    std::cout << tabs << "Workloads:\n    " << tabs;
    for (auto w: workloads)
        std::cout << std::setprecision(4) << w << " ";
    std::cout << "\n";


    double mean = math::mean(workloads);
    std::cout << tabs << "Mean:           " << std::setw(10) << mean << "\n";


    double sigma = math::std_deviation(workloads, mean);
    std::cout << tabs << "Std. deviation: " << std::setw(10) << sigma << " ( " << std::setprecision(2)
              << 100.0 * sigma / mean << " % )\n";

    double diff = 0.0;
    for (auto w: workloads)
        diff = std::max(diff, fabs(w - mean));
    std::cout << tabs << "Max. deviation: " << std::setw(10) << diff << " ( " << std::setprecision(2)
              << 100.0 * diff / mean << " % )\n";
}