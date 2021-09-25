//
// Created by 159-mrv on 6/28/18.
//

#include <utils/math.h>

size_t math::int_log(double x, double base) {
    return static_cast<size_t>(std::log(x)/std::log(base));
}

double math::mean(const vector<double>& values) {
    double avg = 0.0;
    for(auto& d: values) {
        avg += d;
    }
    return avg / values.size();
}

double math::std_deviation(const vector<double>& values, double mean) {
    if (std::isnan(mean)) {
        mean = math::mean(values);
    }
    double sigma = 0.0;
    for (auto & d: values)
        sigma += (d - mean) * (d - mean);

    return sqrt(sigma / (values.size() - 1));
}