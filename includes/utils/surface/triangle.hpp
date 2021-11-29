#pragma once

#include <vector>
#include <string>
#include <map>
#include <set>
#include <utility>
#include <utils/surface/vertex.hpp>
#include <boost/container/static_vector.hpp>
#include <Dense> /// from Eigen library

using boost::container::static_vector;
using Eigen::Vector3d;

namespace surf
{
    /// @class Описывает треугольник
    class Triangle
    {
    public:
        size_t index;

        Triangle() = default;

        explicit Triangle(const std::vector<std::pair<Vertex, size_t>> &verts);

        Triangle(size_t index, const std::vector<size_t> &verts, const Vector3d &n, double area);

        ~Triangle() = default;

        void set_triangles(const std::vector<size_t> &t);

        /// @brief Индексы вершин треугольника
        std::array<size_t, 3> vertices;

        /// @brief Индексы смежных треугольников
        /// смежные треугольники - треугольники имеющие 2 общие вершины с данным
        std::array<size_t, 3> triangles;

        /// @brief Внешняя нормаль треугольника
        Vector3d normal;

        /// @brief площадь треугольника
        double area;
    };
}