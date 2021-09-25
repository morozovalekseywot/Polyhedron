//
// Created by 159-egi on 8/22/17.
//

#ifndef AMR_VERTEX_H
#define AMR_VERTEX_H

#include <allstd.h>

#include <core/indexed_object.h>


class Vertex : public IndexedObject {
public:
    using Ptr = shared_ptr<Vertex>;
    using Ref = const shared_ptr<Vertex>&;

    Vertex(double x, double y, double z = 0.0)
        : IndexedObject(), m_v({x, y, z}) {
    }

    explicit Vertex(const Vector2d& v)
        : IndexedObject(), m_v({v[0], v[1], 0.0}) {
    }

    explicit Vertex(const Vector3d& v)
        : IndexedObject(), m_v(v) {
    }

    inline double x() const noexcept { return m_v[0]; }

    inline double y() const noexcept { return m_v[1]; }

    inline double z() const noexcept { return m_v[2]; }

    inline double r() const noexcept { return m_v.norm(); }

    inline double phi() const noexcept { return std::atan2(m_v[1], m_v[0]); }

    inline void move(const Vector3d& v) { m_v = v; }

    // 2D ATAVISM. Need to remove.
    inline void move_2d(const Vector2d &v) { m_v[0] = v[0]; m_v[1] = v[1]; }

    inline const Vector3d& v() const noexcept { return m_v; }

    // 2D ATAVISM. Need to remove.
    inline Vector2d v_2d() const noexcept { return { m_v[0], m_v[1] }; }

    static inline Vertex::Ptr create(double x, double y, double z = 0.0) {
        return make_shared<Vertex>(x, y, z);
    }

    static inline Vertex::Ptr create_2d(const Vector2d &v) {
        return make_shared<Vertex>(v);
    }

    static inline Vertex::Ptr create(const Vector3d &v) {
        return make_shared<Vertex>(v);
    }

    static inline double distance(Vertex::Ref v1, Vertex::Ref v2) {
        return (v1->m_v - v2->m_v).norm();
    }

    static inline double distance2(Vertex::Ref v1, Vertex::Ref v2) {
        return fabs(v1->m_v[0] - v2->m_v[0]) + fabs(v1->m_v[1] - v2->m_v[1]) + fabs(v1->m_v[2] - v2->m_v[2]);
    }

    inline bool equal(Vertex::Ref v) const {
        return (fabs(m_v[0] - v->m_v[0]) + fabs(m_v[1] - v->m_v[1]) + fabs(m_v[2] - v->m_v[2])) < 1.0e-11;
    }

private:
    Vector3d m_v;
};

/// @brief Структура вершины для сериализации
struct ClearVertex {
    IndexedObject::index_type index;
    double x;
    double y;
    double z;

    ClearVertex() = default;

    explicit ClearVertex(const Vertex::Ptr& vertex) {
        index = vertex->get_index();

        x = vertex->x();
        y = vertex->y();
        z = vertex->z();
    }

    explicit ClearVertex(const string& line) {
        istringstream ss(line);
        ss >> index >> x >> y >> z;
    }
};

/// @brief Вывод ClearVertex в поток в ASCII формате
ostream& operator<<(ostream& out, const ClearVertex& vertex);

#endif //AMR_VERTEX_H
