#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/side.h>
#include <core/cell/cell.h>
#include <core/generator/geometry_generator.h>
#include <core/mesh/mesh.h>
#include <core/cell/data_holder.h>
#include <utils/memory/list_position.h>
#include <utils/geom.h>


// ============================================================================
//                             КОНСТРУКТОРЫ
// ============================================================================

Cell::Cell(Face::Ref left, Face::Ref bottom, Face::Ref right, Face::Ref top)
    : Node(), IndexedObject() {

    m_faces.resize(4);
    m_faces[to_int(Side::LEFT)]   = FacesList::create(left);
    m_faces[to_int(Side::RIGHT)]  = FacesList::create(right);
    m_faces[to_int(Side::BOTTOM)] = FacesList::create(bottom);
    m_faces[to_int(Side::TOP)]    = FacesList::create(top);

    initialize();
}

Cell::Cell(FacesList::Ref left , FacesList::Ref bottom,
           FacesList::Ref right, FacesList::Ref top)
        : Node(), IndexedObject() {

    m_faces.resize(4);
    m_faces[to_int(Side::LEFT)]   = left;
    m_faces[to_int(Side::RIGHT)]  = right;
    m_faces[to_int(Side::BOTTOM)] = bottom;
    m_faces[to_int(Side::TOP)]    = top;

    initialize();
}

Cell::Cell(Face::Ref left, Face::Ref bottom, Face::Ref right,
           Face::Ref top, Face::Ref front, Face::Ref back) {

    m_faces.resize(6);
    m_faces[to_int(Side::LEFT)]   = FacesList::create(left);
    m_faces[to_int(Side::RIGHT)]  = FacesList::create(right);
    m_faces[to_int(Side::BOTTOM)] = FacesList::create(bottom);
    m_faces[to_int(Side::TOP)]    = FacesList::create(top);
    m_faces[to_int(Side::FRONT)]  = FacesList::create(front);
    m_faces[to_int(Side::BACK)]   = FacesList::create(back);

    initialize();
}

Cell::Cell(FacesList::Ref left, FacesList::Ref bottom, FacesList::Ref right,
           FacesList::Ref top, FacesList::Ref front, FacesList::Ref back) {

    m_faces.resize(6);
    m_faces[to_int(Side::LEFT)]   = left;
    m_faces[to_int(Side::RIGHT)]  = right;
    m_faces[to_int(Side::BOTTOM)] = bottom;
    m_faces[to_int(Side::TOP)]    = top;
    m_faces[to_int(Side::FRONT)]  = front;
    m_faces[to_int(Side::BACK)]   = back;

    initialize();
}

void Cell::initialize() {
    m_level = 0;
    m_z = 0;
    m_id = static_cast<IDType>(-1);

    //m_parent = nullptr;
    m_children.clear();

    m_owner = mpi::rank();

    m_data = DataHolder::create(0);

    m_adaptation_flag = AdaptationFlag::NONE;

    m_workload = 0.0;


    // Инициализация геометрии
    Vector3d cntr = {0.0, 0.0, 0.0};
    for (auto &face: m_faces) {
        cntr += face->center();
    }
    m_center = cntr / m_faces.size();

    m_volume = 0.0;
    for (auto &face: m_faces) {
        Vector3d fc = face->center();
        double h = fabs(face->outward_normal().dot(m_center - fc));
        m_volume += h * face->area();
    }
    _2D(m_volume /= 2.0;)
    _3D(m_volume /= 3.0;)

    // Ориентация граней
    size_t n = m_faces.size();
    m_or.resize(n);

    for (size_t i = 0; i < n; ++i) {
        m_or[i] = (m_faces[i]->center() - m_center).dot(m_faces[i]->outward_normal()) > 0.0;
    }

    m_mass_center = {0.0, 0.0, 0.0};
#if DIM2
    auto vertices = ordered_vertices();

    auto j = vertices.size() - 1;

    for (uint i = 0; i < vertices.size(); ++i) {
        double x1 = vertices[j]->x();
        double y1 = vertices[j]->y();

        double x2 = vertices[i]->x();
        double y2 = vertices[i]->y();

        m_mass_center[0] += (y2 - y1) * (x2 * x2 + x2 * x1 + x1 * x1);
        m_mass_center[1] -= (x2 - x1) * (y2 * y2 + y2 * y1 + y1 * y1);

        j = i;
    }

    m_mass_center /= (6.0 * m_volume);
#else
    m_mass_center = m_center;
#endif
}


// ============================================================================
//                              ПРИМИТИВЫ
// ============================================================================

bool Cell::is_correct_face(Side side) const {
    return m_or[to_int(side)];
}

bool Cell::is_correct_face(Face::Ref face) const {
    return face->outward_normal().dot(face->center() - m_center) > 0.0;
}

bool Cell::is_correct_face(FacesList::Ref faces) const {
    return faces->outward_normal().dot(faces->center() - m_center) > 0.0;
}

Side Cell::opposite_side(Side side) const {
    auto faces = m_faces[to_int(side)];
    Face::Ptr face = nullptr;
    for(size_t i = 0; i < faces->size(); ++i) {
        if (faces->at(i)->neighbor_exist(this)) {
            face = faces->at(i);
        }
    }
    if (!face) {
        throw runtime_error("Cell::opposite_side(...) error: cell has no neighbor");
    }

    auto neighbor = face->neighbor(this);

#if DIM2
    auto side_list = TWO_DIMENSIONAL_SIDES;
#else
    auto side_list = THREE_DIMENSIONAL_SIDES;
#endif

    double L = Vertex::distance2(face->vertex(0), face->vertex(1));
    for(auto opp_side: side_list) {
        auto opp_faces = neighbor->faces(opp_side);
        for(size_t i = 0; i < opp_faces->size(); ++i) {
            auto opp_face = neighbor->faces(opp_side)->at(i);

            Vector3d diff = face->center() - opp_face->center();
            if (fabs(diff[0]) + fabs(diff[1]) + fabs(diff[2]) < 1.0e-6 * L) {
                return opp_side;
            }
        }
    }

    std::cerr << "==================== Can't find opposite side ====================\n";
    std::cerr << "Rank: " << mpi::rank() << "\n";
    std::cerr << "Cell.z: " << m_z << ". Side: " << to_string(side) << "\n";
    std::cerr << (is_remote() ? "remote" : "local") << " cell\n";
    std::cerr << "Faces: ";
    for(size_t i = 0; i < faces->size(); ++i) {
        std::cerr << to_string(faces->at(i));
    }
    std::cerr << "Nei.z: " << neighbor->z() << "\n";
    std::cerr << "Nei.faces:\n";
    for(auto& f: neighbor->faces_list()) {
        std::cerr << to_string(f) << "\n";
    }
    std::cerr << "Nei roles:\n";
    std::cerr << (neighbor->get_lists().empty() ? "EMPTY BLYAD" : "NOT EMPTY") << "\n";
    for(auto list: neighbor->get_lists()) {
        std::cerr << "role " << static_cast<int>(list->role()) << "\n";
    }

    throw runtime_error("Cell::opposite_side(...) error: Can't find opposite side");
}

static_vector<Face::Ptr, 24> Cell::faces_list() const {
    static_vector<Face::Ptr, 24> all_faces;
    for (auto& faces: m_faces) {
        for (size_t j = 0; j < faces->size(); ++j) {
            all_faces.emplace_back(faces->at(j));
        }
    }
    return all_faces;
}

static_vector<Face::Ptr, 24> Cell::faces() const {
    static_vector<Face::Ptr, 24> all_faces;
    for (auto& faces: m_faces) {
        for (size_t j = 0; j < faces->size(); ++j) {
            all_faces.emplace_back(faces->at(j));
        }
    }
    return all_faces;
}

FacesList::Ref Cell::faces(Side side) const {
    return m_faces[to_int(side)];
}

void Cell::set_faces(Side side, FacesList::Ref faces) {
    m_center += (faces->center() - m_faces[to_int(side)]->center()) / m_faces.size();
    m_faces[to_int(side)] = faces;
    m_or[to_int(side)] = (faces->center() - m_center).dot(faces->outward_normal()) > 0.0;
}

Vertex::Ref Cell::corner(Side side1, Side side2) const {
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (m_faces[to_int(side1)]->vertex(i)->equal(m_faces[to_int(side2)]->vertex(j))) {
                return m_faces[to_int(side1)]->vertex(i);
            }
        }
    }
    throw runtime_error("Cell::corner(...) error: Can't find corner");
}

Vertex::Ref Cell::left_top_vertex() const noexcept {
    return corner(Side::TOP, Side::LEFT);
}

Vertex::Ref Cell::left_bottom_vertex() const noexcept {
    return corner(Side::BOTTOM, Side::LEFT);
}

Vertex::Ref Cell::right_bottom_vertex() const noexcept {
    return corner(Side::RIGHT, Side::BOTTOM);
}

Vertex::Ref Cell::right_top_vertex() const noexcept {
    return corner(Side::RIGHT, Side::TOP);
}

static_vector<Vertex *, 8> Cell::ordered_vertices() const {
    static_vector<Vertex *, 8> vertices;

    for(auto side : {Side::BOTTOM, Side::RIGHT, Side::TOP, Side::LEFT}) {
        int idx = to_int(side);
        vertices.push_back(m_or[idx] ?
                           m_faces[idx]->vertex_raw(0) :
                           m_faces[idx]->vertex_raw(1));

        if (m_faces[idx]->is_complex()) {
            vertices.push_back(m_faces[idx]->at(0)->vertex_raw(1));
        }
    }

    return vertices;
}

static_vector<Cell*, 24> Cell::neighbors_raw() const {
    static_vector<Cell*, 24> neighbors;
    for(auto& faces: m_faces) {
        faces->append_neighbors(this, neighbors);
    }
    return neighbors;
}

static_vector<Cell::Ptr, 24> Cell::neighbors() const {
    static_vector<Cell::Ptr, 24> neighbors;
    for(auto& faces: m_faces) {
        faces->append_neighbors(this, neighbors);
    }
    return neighbors;
}


// ============================================================================
//                              ИНДЕКСЫ
// ============================================================================

bool Cell::equal_to(Cell::Ref cell) const noexcept {
    return cell->level() == level() && cell->m_z == m_z;
}

size_t Cell::level() const noexcept {
    return m_level;
}

void Cell::set_level(size_t lvl) {
    m_level = lvl;
}

size_t Cell::z_base() const noexcept {
    return m_z / math::pow4(m_level);
}

size_t Cell::z() const noexcept {
    return m_z;
}

void Cell::set_z(size_t z) {
    m_z = z;
}

Cell::IDType Cell::id() const noexcept {
    return m_id;
}

void Cell::set_id(IDType id) {
    m_id = id;
}


// ============================================================================
//                              ИЕРАРХИЯ
// ============================================================================

Cell::Ptr Cell::parent() const noexcept {
    return m_parent.lock();
}

void Cell::set_parent(Cell::Ref parent) {
    m_parent = parent;
}

Cell::Ptr Cell::origin(Cell::Ptr cell)  {
    return cell->m_parent.lock() ? Cell::origin(cell->m_parent.lock()) : cell;
}

const vector<Cell::Ptr>& Cell::children() const noexcept {
    return m_children;
}

void Cell::leaf_children(Cell::Ref cell, vector<Cell::Ptr>& leaf) {
    if (cell->m_children.empty()) {
        leaf.emplace_back(cell);
    }
    else {
        for(auto& child: cell->m_children) {
            Cell::leaf_children(child, leaf);
        }
    }
}

static_vector<Cell::Ptr, 4> Cell::children_by_side(Side side) const {
#if DIM2
    switch (side) {
        case Side::LEFT:
            return {m_children[0], m_children[2]};
        case Side::BOTTOM:
            return {m_children[0], m_children[1]};
        case Side::RIGHT:
            return {m_children[1], m_children[3]};
        case Side::TOP:
            return {m_children[2], m_children[3]};
        default:
            throw runtime_error("Error! Rectangular cell doesn't contain that size");
    }
#else
    switch(side) {
        case Side::LEFT:
            return { m_children[0], m_children[1], m_children[3], m_children[2] };
        case Side::RIGHT:
            return { m_children[4], m_children[5], m_children[7], m_children[6] };
        case Side::BOTTOM:
            return { m_children[0], m_children[1], m_children[5], m_children[4] };
        case Side::TOP:
            return { m_children[2], m_children[3], m_children[7], m_children[6] };
        case Side::BACK:
            return { m_children[0], m_children[2], m_children[6], m_children[4] };
        case Side::FRONT:
            return { m_children[1], m_children[3], m_children[7], m_children[5] };
        default:
            throw runtime_error("Error! Cubic cell doesn't contain that size");
    }
#endif
}

Cell::Ref Cell::child_directed_to(Cell::Ref cell) const {
    auto level_dif = cell->level() - m_level;
    size_t dz = cell->z() / math::pow4(level_dif - 1);
    size_t loc_z = dz % 4;
    return get_child_by_local_id(loc_z);
}

void Cell::set_children(const vector<Cell::Ptr>& children) {
    m_children = children;
}

Cell::Ref Cell::get_child_by_local_id(size_t loc_z) const {
    if (loc_z > m_children.size()) {
        throw runtime_error("Error: Trying to get child with local ID more than number of children");
    }

    return m_children[loc_z];
}

void Cell::set_child_by_local_id(Cell::Ref cell, size_t loc_z, Cell::Ref child) {
    if (loc_z > cell->m_children.size()) {
        throw runtime_error("Error: Trying to get child with local ID more than number of children");
    }

    cell->m_children[loc_z] = child;
    child->set_parent(cell);
}

bool Cell::is_childless() const noexcept {
    return m_children.empty();
}

void Cell::remove_children() {
    for (auto &ch: m_children) {
        ch->remove_children();
    }
    m_children.clear();
}

static_vector<Cell*, 7> Cell::siblings() const {
    static_vector<Cell*, 7> sibs;
    if (!m_parent.lock()) {
        return sibs;
    }
    else {
        for(Cell::Ref child: m_parent.lock()->m_children) {
            if (child.get() != this) {
                sibs.push_back(child.get());
            }
        }
        return sibs;
    }
}


// ============================================================================
//                              ДАННЫЕ, ФЛАГИ
// ============================================================================

int Cell::owner() const noexcept {
    return m_owner;
}

void Cell::set_owner(int rank) {
    m_owner = rank;
    if (m_parent.lock()) {
        m_parent.lock()->set_owner(rank);
    }
}

bool Cell::is_local() const noexcept {
    return m_owner == mpi::rank();
}

bool Cell::is_remote() const noexcept {
    return m_owner != mpi::rank();
}

DataHolder::Ref Cell::data_holder() const {
    return m_data;
}

double * Cell::data() const {
    return m_data->buffer();
}

void Cell::set_data(DataHolder::Ref data) {
    m_data = data;
}

void Cell::copy_data_from(Cell::Ref cell) {
    m_data->resize(cell->m_data->size());
    m_data->deserialize(cell->m_data->buffer());
}

AdaptationFlag Cell::adaptation_flag() const noexcept {
    return m_adaptation_flag;
}

void Cell::set_adaptation_flag(AdaptationFlag flag) {
    m_adaptation_flag = flag;
}

bool Cell::can_coarse() const {
    // Установим детей
    if (m_children.empty()) {
        return m_adaptation_flag == AdaptationFlag::COARSE;
    }

    // Если ячейка локальная, тогда все дети должны быть готовы к огрублению
    if (is_local()) {
        for(auto& child: m_children) {
            if (!child->can_coarse()) {
                return false;
            }
        }
        return true;
    }
    else {
        // Если ячейка удаленная, достаточно двух детей готовых к огрублению
        int children_to_coarse = 0;

        for(auto& child: m_children) {
            if (child->can_coarse()) {
                ++children_to_coarse;
            }
        }
        return children_to_coarse > 1;
    }
}

double Cell::workload() const noexcept {
    return m_workload;
}

void Cell::reset_workload(){
    m_workload = 0.0;
}

void Cell::add_workload(double workload) {
    m_workload += workload;
}


// ============================================================================
//                                 ГЕОМЕТРИЯ
// ============================================================================


bool Cell::is_2D() const noexcept {
    return m_faces.size() < 5;
}

bool Cell::is_3D() const noexcept {
    return m_faces.size() > 5;
}

const Vector3d & Cell::center() const noexcept {
    return m_center;
}

const Vector3d & Cell::mass_center() const noexcept {
    return m_mass_center;
}

Vector2d Cell::center_2d() const noexcept {
    return { m_center[0], m_center[1] };
}

double Cell::center_radius() const {
    return m_center.norm();
}

double Cell::center_angle() const {
    return atan2(m_center[1], m_center[0]);
}

double Cell::volume() const {
    return m_volume;
}

double Cell::volume_as() const {
    double G = 0.0;

    auto vertices = ordered_vertices();

    auto j = vertices.size() - 1;

    for (uint i = 0; i < vertices.size(); ++i) {
        double x1 = vertices[j]->x();
        double y1 = vertices[j]->y();

        double x2 = vertices[i]->x();
        double y2 = vertices[i]->y();

        G += (x2 - x1) * (y2 * y2 + y2 * y1 + y1 * y1);

        j = i;
    }

    G = -G / 6.0;

    if (G < 0) {
        std::cerr << "Error! Neg volume: " << G << "\n";
        throw runtime_error("ERROR! Neg volume");
    }

    return G;
}

Vector3d Cell::volume_as(const Vector3d & v0) const {
    double x0(v0[0]), y0(v0[1]);
    double vol_x(0.0), vol_y(0.0), vol_z(0.0);

    auto vertices = ordered_vertices();

    auto j = vertices.size() - 1;

    for (uint i = 0; i < vertices.size(); ++i) {
        double x1 = vertices[j]->x() - x0;
        double y1 = vertices[j]->y() - y0;

        double x2 = vertices[i]->x() - x0;
        double y2 = vertices[i]->y() - y0;

        vol_x += (y2 - y1) * (x2 * x2 + x2 * x1 + x1 * x1);
        vol_y += (x2 - x1) * (y2 * y2 + y2 * y1 + y1 * y1);

        j = i;
    }

    vol_x = vol_x / 6.0;
    vol_y = -vol_y / 6.0;

    return {vol_x, vol_y, vol_z};
}

// ============================================================================
//                                 ДЕБАГ
// ============================================================================

void Cell::print_dbg_info(bool neighbors_out) const {
    std::cout << "Cell:\n (Z,LVL) => (" << "," << level() << ")\n Faces: ";

    std::cout << " Remote: " << is_remote() << std::endl;

    if(neighbors_out) {
        std::cout << "; \n Neighbors: \n";

        for(auto& neighbor: neighbors_raw()) {
            neighbor->print_dbg_info(false);
        }

        std::cout << std::endl;
    }

    std::cout << " Data: " << "\n";
    for(size_t q = 0; q < 5; ++q) {
        std::cout << "  quant " << q << ": ";
        for(uint i = 0; i < m_data->size()/5; ++i) {
            std::cout << m_data->buffer()[5 * q + i] << " ";
        }
        std::cout << "\n";
    }
}

bool Cell::interesting() const {
    if (m_z == 3200 || m_z == 3299) {
        return true;
    }
    else {
        return false;
    }
}

ClearCell::ClearCell(Cell::Ref cell, bool is_fake) {
    index = cell->get_index();

    id = cell->id();
    z = cell->z();
    level = cell->level();

    owner = cell->owner();
    fake = is_fake;

    for (int s = 0; s < 4; ++s) {
        auto side = to_side(s);

        auto lol_faces = cell->faces(side);
        faces[2 * s] = lol_faces->at(0)->get_index();
        faces[2 * s + 1] = lol_faces->is_complex() ?
                           lol_faces->at(1)->get_index() :
                           IndexedObject::undefined_index;
    }
}

ClearCell::ClearCell(const string &line) {
    std::istringstream ss(line);

    ss >> index >> id >> z >> level >> owner;
    ss >> faces[0] >> faces[1] >> faces[2] >> faces[3];
    ss >> faces[4] >> faces[5] >> faces[6] >> faces[7];
    ss >> fake;
}

ostream& operator<<(ostream& out, const ClearCell& cell) {
    out << cell.index << " " << cell.id << " " << cell.z << " "
        << cell.level << " " << cell.owner << " "
        << cell.faces[0] << " " << cell.faces[1] << " " << cell.faces[2] << " " << cell.faces[3] << " "
        << cell.faces[4] << " " << cell.faces[5] << " " << cell.faces[6] << " " << cell.faces[7] << " "
        << cell.fake;
    return out;
}
