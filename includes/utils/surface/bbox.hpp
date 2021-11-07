#pragma once

#include <utils/surface/face.hpp>

// TODO разделить по файлам
namespace surf
{
    struct BBox
    {
        Eigen::Vector3d x0, len;

        explicit BBox(const std::vector<Face> &faces)
        {
            double xmin = 1e5, xmax = -1e5, ymin = 1e5, ymax = -1e5, zmin = 1e5, zmax = -1e5;
            for (auto &face: faces)
            {
                for (auto &v: face.vects)
                {
                    xmin = min(xmin, v.x());
                    xmax = max(xmax, v.x());
                    ymin = min(ymin, v.y());
                    ymax = max(ymax, v.y());
                    zmin = min(zmin, v.z());
                    zmax = max(zmax, v.z());
                }
            }

            x0 = Eigen::Vector3d(xmin, ymin, zmin);
            len = Eigen::Vector3d(xmax - xmin, ymax - ymin, zmax - zmin);

            // расширение BBox
//            x0 = x0 + Eigen::Vector3d(-abs(xmin), -abs(ymin), -abs(zmin)) * 0.07;
//            // xmax + abs(xmax)*0.02 - (xmin - abs(xmin)*0.02)
//            len = len + Eigen::Vector3d(abs(xmax) + abs(xmin), abs(ymax) + abs(ymin), abs(zmax) + abs(zmin)) * 0.07;

        }

        BBox(const BBox &bbox)
        {
            x0 = bbox.x0;
            len = bbox.len;
        }

        double volume() const
        {
            return len.x() * len.y() * len.z();
        }

        ~BBox() = default;
    };
}

