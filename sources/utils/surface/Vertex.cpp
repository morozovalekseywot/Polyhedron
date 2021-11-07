#include <utils/surface/Vertex.hpp>

surf::Vertex::Vertex(const Vector3d &point) : v(point)
{}

Vector3d surf::Vertex::operator+(const surf::Vertex &a) const
{
    return a.v + v;
}

Vector3d surf::Vertex::operator-(const surf::Vertex &a) const
{
    return v - a.v;
}

double surf::Vertex::mod() const
{
    return v.norm();
}

void surf::Vertex::normalize()
{
    v.normalize();
}

template<class T>
Vector3d surf::Vertex::operator*(const T &a) const
{
    return v * a;
}

template<class T>
Vector3d surf::Vertex::operator/(const T &a) const
{
    return v / a;
}

double surf::Vertex::operator*(const surf::Vertex &a) const
{
    return a.v.x() * v.x() + a.v.y() * v.y() + a.v.z() * v.z();
}

Vector3d surf::Vertex::operator-()
{
    return -v;
}

void surf::Vertex::set_next_triangle(size_t idx)
{
    if (triangles.size() >= max_triangles)
        throw std::runtime_error("number_of_triangles > max_triangles");
    triangles.push_back(idx);
}

Vector3d cross(const surf::Vertex &a, const surf::Vertex &b)
{
    return a.v.cross(b.v);
}

std::ostream &operator<<(std::ostream &os, const surf::Vertex &a)
{
    os << a.v.x() << " " << a.v.y() << " " << a.v.z();
    return os;
}
