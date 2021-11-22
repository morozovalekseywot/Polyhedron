#pragma once

#include <utils/surface/face.hpp>

namespace surf
{
    struct BBox
    {
        Eigen::Vector3d x0, len;

        explicit BBox(const std::vector<Face> &faces);

        BBox(const BBox &bbox);

        double volume() const;

        ~BBox() = default;
    };
}

