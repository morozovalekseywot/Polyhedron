//
// Created by 159-mrv on 5/21/18.
//

#ifndef INFINIMESH_ALLSTD_H
#define INFINIMESH_ALLSTD_H

// Define
#define DIM DIMENSION
#if DIM == 2
    #define DIM2 1
    #define DIM3 0
    #define _2D(...) __VA_ARGS__
    #define _3D(...)
#else
    #define DIM2 0
    #define DIM3 1
    #define _2D(...)
    #define _3D(...) __VA_ARGS__
#endif



// STD
#include <memory>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <utility>
#include <atomic>
#include <thread>
#include <mutex>
#include <limits>
#include <cmath>
#include <stdexcept>

using std::shared_ptr;
using std::weak_ptr;
using std::unique_ptr;
using std::make_shared;
using std::runtime_error;
using std::min;
using std::max;
using std::swap;
using std::isnan;
using std::thread;
using std::get;
using std::tie;

// define makeUnique
// using std::make_unique;
template<typename T, typename... Args>
std::unique_ptr<T> makeUnique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
};

#ifndef __USE_MISC
typedef unsigned int uint;
#endif

// STD containers
#include <vector>
#include <string>
#include <array>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>

using std::vector;
using std::string;
using std::array;
using std::unordered_set;
using std::unordered_map;
using std::tuple;

using std::make_pair;
using std::to_string;


// STD streams
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

using std::ifstream;
using std::ofstream;
using std::ostream;
using std::istream;
using std::stringstream;
using std::ostringstream;
using std::istringstream;

using std::setprecision;
using std::scientific;
using std::fixed;


// Boost
#include <boost/container/static_vector.hpp>

using boost::container::static_vector;


// LIBS
#include <Dense>

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix3d;
using Eigen::Matrix2d;


#endif //INFINIMESH_ALLSTD_H
