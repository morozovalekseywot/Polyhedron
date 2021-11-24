#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <boost/container/static_vector.hpp>
#include <Dense> /// from Eigen library

namespace surf
{
    struct Face
    {
        const size_t n_verts = 3;

        std::array<Eigen::Vector3d, 3> vects; /// массив вершин
        Eigen::Vector3d normal;        /// единичная внешняя нормаль

        Face(const std::array<Eigen::Vector3d,3> &a, const Eigen::Vector3d &norm);

        explicit Face(const std::vector<Eigen::Vector3d> &a);

        Face &operator=(const Face &a);
    };
}