//
// Created by 159-mrv on 8/22/19.
//

#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/cell/cell.h>
#include <utils/geom.h>

Face::Face(Vertex::Ref begin, Vertex::Ref end)
        : m_vertices({begin, end}),
          m_flag(FaceFlag::ORDER) {
    initialize();
}

#if DIM2
void Face::initialize() {
    if (m_vertices[0] && m_vertices[1]) {
        m_center = 0.5 * (m_vertices[0]->v() + m_vertices[1]->v());

        m_area   = (m_vertices[1]->v() - m_vertices[0]->v()).norm();

        Vector3d face_vec = (m_vertices[1]->v() - m_vertices[0]->v()) / m_area;

        m_outward_normal = { face_vec[1], -face_vec[0], 0.0 };
        m_inward_normal  = -m_outward_normal;
    } else {
        throw runtime_error("Face::initialize_1d() error: Attempt to create face on nullptr vertices");
    }
}
#endif

Face::Face(Vertex::Ref v1, Vertex::Ref v2, Vertex::Ref v3, Vertex::Ref v4)
        : m_vertices({v1, v2, v3, v4}),
          m_flag(FaceFlag::ORDER) {
    initialize();
}

#if DIM3
void Face::initialize() {
    if (m_vertices[0] && m_vertices[1] && m_vertices[2] && m_vertices[3]) {
        array<Vector3d, 4> vs = {
            m_vertices[0]->v(),
            m_vertices[1]->v(),
            m_vertices[2]->v(),
            m_vertices[3]->v()
        };

        m_center = 0.25 * (vs[0] + vs[1] + vs[2] + vs[3]);

        m_area = geom::polygon_area(vs);

        m_outward_normal = geom::normal(vs);
        m_inward_normal  = -m_outward_normal;
    } else {
        throw runtime_error("Face::initialize_2d() error: Attempt to create face on nullptr vertices");
    }
}
#endif

void Face::reverse() {
#if DIM2
        std::swap(m_vertices[0], m_vertices[1]);
#else
        std::swap(m_vertices[0], m_vertices[3]);
        std::swap(m_vertices[1], m_vertices[2]);
#endif
    std::swap(m_inward_normal, m_outward_normal);
    std::swap(m_inward_cell, m_outward_cell);
}

#if DIM3
void Face::rotate(int r) {
    auto v0 = m_vertices[(4 - r) % 4];
    auto v1 = m_vertices[(5 - r) % 4];
    auto v2 = m_vertices[(6 - r) % 4];
    auto v3 = m_vertices[(7 - r) % 4];
    m_vertices = {v0, v1, v2, v3};
}
#endif

double Face::area() const {
    return m_area;
}

double Face::area_as() const {
    return m_area * m_center[1];
}

const Vector3d& Face::center() const {
    return m_center;
}

Vector2d Face::center_2d() const {
    return {m_center[0], m_center[1] };
}

const Vector3d& Face::normal(Cell::Ref neighbor) const {
    if(neighbor == m_inward_cell.lock()) {
        return m_inward_normal;
    }
    else if(neighbor == m_outward_cell.lock()) {
        return m_outward_normal;
    }
    else {
        throw runtime_error("Face::normal(...) error: Attempt to get normal to not adjacent neighbor");
    }
}

const Vector3d& Face::normal(Cell* neighbor) const {
    if(neighbor == m_inward_cell.lock().get()) {
        return m_inward_normal;
    }
    else if(neighbor == m_outward_cell.lock().get()) {
        return m_outward_normal;
    }
    else {
        throw runtime_error("Face::normal(...) error: Attempt to get normal to not adjacent neighbor");
    }
}

Vector2d Face::normal_2d(Cell::Ref neighbor) const {
    Vector3d n = normal(neighbor);
    return { n[0], n[1] };
}

const Vector3d & Face::outward_normal() const {
    return m_outward_normal;
}

FaceFlag Face::flag() const {
    return m_flag;
}

void Face::set_flag(FaceFlag flag) {
    m_flag = flag;
}

const vector<Vertex::Ptr>& Face::vertices() const noexcept {
    return m_vertices;
}

size_t Face::size() const noexcept {
    return m_vertices.size();
}

Vertex::Ref Face::vertex(int idx) const noexcept {
    return m_vertices[idx];
}
void Face::set_vertex(int idx, Vertex::Ref v) noexcept {
    m_vertices[idx] = v;
    initialize();
}

Vertex* Face::vertex_raw(int idx) const noexcept {
    return m_vertices[idx].get();
}

Cell::Ptr Face::inward_cell() const {
    return m_inward_cell.lock();
}

void Face::set_inward_cell(Cell::Ref cell) {
    m_inward_cell = cell;
}

Cell::Ptr Face::outward_cell() const {
    return m_outward_cell.lock();
}

void Face::set_outward_cell(Cell::Ref cell) {
    m_outward_cell = cell;
}

bool Face::neighbor_exist(const Cell * cell) const noexcept {
    if (cell == m_inward_cell.lock().get()) {
        if (m_outward_cell.lock()) {
            return true;
        } else {
            return false;
        }
    } else if (cell == m_outward_cell.lock().get()) {
        if (m_inward_cell.lock()) {
            return true;
        } else {
            return false;
        }
    }
    return false;
}

bool Face::neighbor_exist(Cell::Ref cell) const noexcept {
    if(cell == m_inward_cell.lock()) {
        if(m_outward_cell.lock()) {
            return true;
        }
        else {
            return false;
        }
    }
    else if(cell == m_outward_cell.lock()) {
        if(m_inward_cell.lock()) {
            return true;
        }
        else {
            return false;
        }
    }
    return false;
}

Cell* Face::neighbor_raw(const Cell * cell) const {
#ifdef FULL_CHECK
    Cell* order_lock   = m_inward_cell.lock().get();
    Cell* opposed_lock = m_outward_cell.lock().get();

    if (cell == order_lock) {
        if (opposed_lock) {
            return opposed_lock;
        } else {
            std::cout << to_string(this) << "\n";
            throw runtime_error("ERROR! Trying to access NULL neighbor (opposed cell)");
        }
    } else if (cell == opposed_lock) {
        if (order_lock) {
            return order_lock;
        } else {
            throw runtime_error("ERROR! Trying to access NULL neighbor (order cell)");
        }
    } else {
        throw runtime_error("ERROR! Trying to access neighbor of NULL cell (order or opposed)");
    }
#else
    return cell == m_inward_cell.lock().get() ?
           m_outward_cell.lock().get() :
           m_inward_cell.lock().get();
#endif
}

Cell::Ptr Face::neighbor(const Cell * cell) const {
#ifdef FULL_CHECK
    if (cell == m_inward_cell.lock().get()) {
        if (m_outward_cell.lock()) {
            return m_outward_cell.lock();
        } else {
            std::cout << "face: " << to_string(this) << "\n";
            std::cout << "cell.center: (" << cell->center_2d()[0] << ", " << cell->center_2d()[1] << ")\n";
            throw runtime_error("ERROR! Trying to access NULL neighbor (opposed cell)");
        }
    } else if (cell == m_outward_cell.lock().get()) {
        if (m_inward_cell.lock()) {
            return m_inward_cell.lock();
        } else {
            throw runtime_error("ERROR! Trying to access NULL neighbor (order cell)");
        }
    } else {
        throw runtime_error("ERROR! Trying to access neighbor of NULL cell (order or opposed)");
    }
#else
    return cell == m_inward_cell.lock().get() ?
           m_outward_cell.lock() :
           m_inward_cell.lock();
#endif
}

Cell::Ptr Face::neighbor(Cell::Ref cell) const {
#ifdef FULL_CHECK
    if(cell == m_inward_cell.lock()) {
        if(m_outward_cell.lock()) {
            return m_outward_cell.lock();
        }
        else {
            std::cout << to_string(this) << "\n";
            throw runtime_error("ERROR! Trying to access NULL neighbor (opposed cell)");
        }
    }
    else if(cell == m_outward_cell.lock()) {
        if(m_inward_cell.lock()) {
            return m_inward_cell.lock();
        }
        else {
            throw runtime_error("ERROR! Trying to access NULL neighbor (order cell)");
        }
    }
    else {
        std::cerr << "Bug cell, z: " << cell->z() << "\n";
        throw runtime_error("ERROR! Trying to access neighbor of NULL cell (order or opposed)");
    }
#else
    return cell == m_inward_cell.lock() ?
           m_outward_cell.lock() :
           m_inward_cell.lock();
#endif
}

string to_string(const Face * face) {
    stringstream ss;

#if DIM2
    ss << "(" << face->vertex(0)->x() << ", " << face->vertex(0)->y() << ") - "
       << "(" << face->vertex(1)->x() << ", " << face->vertex(1)->y() << ")";
    return ss.str();

#else
    ss << "(" << face->vertex(0)->x() << ", " << face->vertex(0)->y() << ", " << face->vertex(0)->z() << ") - "
       << "(" << face->vertex(1)->x() << ", " << face->vertex(1)->y() << ", " << face->vertex(1)->z() << ") - "
       << "(" << face->vertex(2)->x() << ", " << face->vertex(2)->y() << ", " << face->vertex(2)->z() << ") - "
       << "(" << face->vertex(3)->x() << ", " << face->vertex(3)->y() << ", " << face->vertex(3)->z() << ")";

    return ss.str();
#endif
}

string to_string(Face::Ref face) {
    return to_string(face.get());
}

ClearFace::ClearFace(Face::Ref face) {
    index = face->get_index();

    vertex[0] = face->vertex(0)->get_index();
    vertex[1] = face->vertex(1)->get_index();
#if DIM2
        vertex[2] = IndexedObject::undefined_index;
        vertex[3] = IndexedObject::undefined_index;
#else
        vertex[2] = face->vertex(2)->get_index();
        vertex[3] = face->vertex(3)->get_index();
#endif

    inward_cell  = face->inward_cell()  ? face->inward_cell()->get_index()  : IndexedObject::undefined_index;
    outward_cell = face->outward_cell() ? face->outward_cell()->get_index() : IndexedObject::undefined_index;

    flag = to_int(face->flag());
}

ClearFace::ClearFace(const string& line) {
    std::istringstream ss(line);
    ss >> index >> vertex[0] >> vertex[1] >> vertex[2] >> vertex[3]
       >> inward_cell >> outward_cell >> flag;
}

ostream& operator<<(ostream& out, const ClearFace& face) {
    out << face.index << " "
        << face.vertex[0] << " " << face.vertex[1] << " "
        << face.vertex[2] << " " << face.vertex[3] << " "
        << face.inward_cell << " " << face.outward_cell << " "
        << face.flag;
    return out;
}