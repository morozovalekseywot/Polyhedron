#include <utils/surface/face.hpp>


surf::Face::Face(const std::array<Eigen::Vector3d,3> &a, const Eigen::Vector3d &norm) : vects(a)
{
    if (a.size() != n_verts)
    {
        throw std::runtime_error("Wrong number of vertices");
    }
    normal = norm;
}

surf::Face::Face(const std::vector<Eigen::Vector3d> &a)
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

surf::Face &surf::Face::operator=(const surf::Face &a)
{
    vects = a.vects;
    normal = a.normal;
    return *this;
}
