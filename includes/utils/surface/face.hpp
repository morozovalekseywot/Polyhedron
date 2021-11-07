#pragma once

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

        std::array<Eigen::Vector3d, 3> vects; // массив вершин
        Eigen::Vector3d normal;        // единичная внешняя нормаль

        Face(const vector<Eigen::Vector3d> &a, const Eigen::Vector3d &norm)
        {
            if (a.size() != n_verts)
            {
                throw std::runtime_error("Wrong number of vertices");
            }
            vects = {a[0], a[1], a[2]};
            normal = norm;
        }

        explicit Face(const vector<Eigen::Vector3d> &a)
        {
            if (a.size() != n_verts)
            {
                throw std::runtime_error("Wrong number of vertices");
            }

            vects = {a[0], a[1], a[2]};

            normal = (vects[1] - vects[0]).cross(vects[2] - vects[0]);
            normal.normalize();
            if (normal == Eigen::Vector3d(0, 0, 0))
            {
                std::cerr << "\n" << vects[0] << "|" << vects[1] << "|" << vects[2] << "\n";
                Eigen::Vector3d vec1 = vects[1] - vects[0];
                Eigen::Vector3d vec2 = vects[2] - vects[0];
                vec1.normalize();
                vec2.normalize();
                std::cerr << vec1 << "|" << vec2;
                throw std::runtime_error("Wrong number of difficult Vertex in face");
            }
        }

        Face &operator=(const Face &a)
        {
            vects = a.vects;
            normal = a.normal;
            return *this;
        }
    };
}