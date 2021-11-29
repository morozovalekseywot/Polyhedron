#include <utils/surface/triangle.hpp>

surf::Triangle::Triangle(const std::vector<std::pair<Vertex, size_t>> &verts)
{
    triangles = {verts[0].second, verts[1].second, verts[2].second};
    normal = (verts[1].first - verts[0].first).cross(verts[2].first - verts[0].first);
    area = 0.5 * ((verts[2].first.v - verts[0].first.v).cross(verts[1].first.v - verts[0].first.v)).norm();
}

surf::Triangle::Triangle(size_t index, const std::vector<size_t> &verts, const Vector3d &n, double area) : index(index), normal(n), area(area)
{
    vertices = {verts[0], verts[1], verts[2]};
}

void surf::Triangle::set_triangles(const std::vector<size_t> &t)
{
    triangles = {t[0], t[1], t[2]};
}
