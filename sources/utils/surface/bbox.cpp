#include <utils/surface/bbox.hpp>

surf::BBox::BBox(const std::vector<Face> &faces)
{
    double xmin = 1e5, xmax = -1e5, ymin = 1e5, ymax = -1e5, zmin = 1e5, zmax = -1e5;
    for (auto &face: faces)
    {
        for (auto &v: face.vects)
        {
            xmin = std::min(xmin, v.x());
            xmax = std::max(xmax, v.x());
            ymin = std::min(ymin, v.y());
            ymax = std::max(ymax, v.y());
            zmin = std::min(zmin, v.z());
            zmax = std::max(zmax, v.z());
        }
    }

    x0 = Eigen::Vector3d(xmin, ymin, zmin);
    len = Eigen::Vector3d(xmax - xmin, ymax - ymin, zmax - zmin);
}

surf::BBox::BBox(const surf::BBox &bbox)
{
    x0 = bbox.x0;
    len = bbox.len;
}

double surf::BBox::volume() const
{
    return len.x() * len.y() * len.z();
}
