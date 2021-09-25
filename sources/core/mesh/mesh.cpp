//
// Created by 159-egi on 8/24/17.
//

#include <control/configuration.h>
#include <control/mpi_wrapper.h>

#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/data_holder.h>
#include <core/layer/border_control.h>
#include <core/layer/ghost_control.h>
#include <core/mesh/mesh.h>
#include <core/generator/sector_generator.h>
#include <core/generator/decomposition.h>
#include <future>

#include <io/vtk/vtk_writer.h>

#include <problems/problem.h>

#include <utils/memory/kd_tree.h>
#include <core/generator/rectangle_generator.h>


// ============================================================================
//                             ИНИЦИАЛИЗАЦИЯ
// ============================================================================

Mesh::Mesh(const Configuration &config) {
    read_configuration(config);
    build(config);
}

unique_ptr<Mesh> Mesh::create(const Configuration &config) {
    return makeUnique<Mesh>(config);
}

Mesh::Mesh(const Configuration &config, const string& checkpoint) {
    read_configuration(config);
    build(config);
    restart_from_checkpoint(checkpoint);
}

unique_ptr<Mesh> Mesh::create(const Configuration &config, const string& checkpoint) {
    return makeUnique<Mesh>(config, checkpoint);
}

Mesh::~Mesh() = default;

void Mesh::read_configuration(const Configuration &config) {
    // Создаем геометрию сетки
    m_geometry = GeometryGenerator::create(config);

    // Врубаем декомпозицию
    m_decomposition = Decomposition::create(config, m_geometry.get());

    // Если существует опция касающаяся адаптации
    // И если опция установлена в положительное значение
    if(config.exist("mesh", "adaptation") && config("mesh", "adaptation").to_bool()) {

        // Устанавливаем максимальный уровень адаптации и пороги
        m_max_level = config("mesh", "adaptation_criteria", "max_fine_level").to_uint();

        // Устанавливаем геометрический тип адаптации
        m_adaptation_type = AdaptationCriterion::NONE;
        if (config.exist("mesh", "adaptation_criteria", "gtype")) {
            string ad_type = config("mesh", "adaptation_criteria", "gtype").to_string();
            if (ad_type == "none") {
                m_adaptation_type = AdaptationCriterion::NONE;
            } else if (ad_type == "all") {
                m_adaptation_type = AdaptationCriterion::ALL;
            } else if (ad_type == "cascade") {
                m_adaptation_type = AdaptationCriterion::CASCADE;
            } else {
                throw runtime_error("Unknown adaptation gtype '" + ad_type + "' in config file.");
            }
        }

        m_solver_adaptation = false;
        if (config.exist("mesh", "adaptation_criteria", "stype")) {
            m_solver_adaptation = config("mesh", "adaptation_criteria", "stype").to_string() != "none";
        }
    } else {
        // Если адаптация не предусмотрена - фиксируем максимальный уровень базовым
        m_max_level = 0;
        m_adaptation_type = AdaptationCriterion::NONE;
        m_solver_adaptation = false;
    }

    m_balancing_options = BalancingOptions(config);

    // Число тредов
    m_n_chunks = 1;
    if (config.exist("calculation", "shared_memory") &&
        config("calculation", "shared_memory").to_bool()) {
        if (config.exist("calculation", "threads_num")) {
            m_n_chunks = config("calculation", "threads_num").to_uint();
        } else {
            m_n_chunks = std::thread::hardware_concurrency();
        }
    }
}

void Mesh::build(const Configuration &config) {
    mpi::cout << "  Start mesh creation\n";

    mpi::cout << "      create and connect base cells\n";
    auto cells = m_geometry->create_cells(m_decomposition.get());

    m_cells = NodeList::create(mpi::rank(), ListRole::MESH);
    for(auto& cell: cells) {
        if (cell->is_local()) {
            m_cells->push_back(cell);
        }
        m_base_cells[cell->z_base()] = cell;
    }

    //dbg_write(nullptr);

    mpi::cout << "      create fake and ghost cells\n";
    auto bc = config("mesh", "boundary_conditions");
    array<FaceFlag, 6> face_flags = {
            to_face_flag(bc(to_string(to_side(0))).to_string()),
            to_face_flag(bc(to_string(to_side(1))).to_string()),
            to_face_flag(bc(to_string(to_side(2))).to_string()),
            to_face_flag(bc(to_string(to_side(3))).to_string()),
            FaceFlag::ORDER,
            FaceFlag::ORDER
    };
    if (m_cells->begin()->is_3D()) {
        if (bc.exist("back") && bc.exist("front")) {
            face_flags[to_int(Side::BACK)] = to_face_flag(bc("back").to_string());
            face_flags[to_int(Side::FRONT)] = to_face_flag(bc("front").to_string());

        } else {
            throw runtime_error("Config does not include 'back' and 'front' boundary conditions");
        }
    }
    initialize_other_cells(face_flags);


    mpi::cout << "      create chunks\n";
    m_cells->shuffle();
    m_cells->split(m_n_chunks);

    for (auto &border: m_border) { border->rebuild(); }
    for (auto &ghost:  m_ghosts) { ghost->rebuild(); }

    for (auto& ghost: m_ghosts) {
        for (auto g_cell: *ghost->cells()) {
            for (auto& face: g_cell->faces_list()) {
                // Удаляем ненужные указатели с гост-ячеек
                if (face->neighbor_exist(g_cell)) {
                    if (face->neighbor(g_cell)->is_remote()) {
                        face->set_inward_cell(nullptr);
                        face->set_outward_cell(nullptr);
                    }
                } else {
                    face->set_inward_cell(nullptr);
                    face->set_outward_cell(nullptr);
                }

                // Восстанавливаем указатель на себя
                if ((face->center() - g_cell->center()).dot(face->outward_normal()) > 0.0) {
                    face->set_inward_cell(g_cell);
                }
                else {
                    face->set_outward_cell(g_cell);
                }
            }
        }
    }

    //geometric_adaptation(m_adaptation_type);

    mpi::cout << "  Base mesh has been created\n" << std::endl;
}

void Mesh::initialize_other_cells(array<FaceFlag, 6> face_flags) {
    m_fake = NodeList::create(mpi::rank());
    std::map<int, NodeList::Ptr> borders;
    std::map<int, NodeList::Ptr> ghosts;

_2D(auto side_list = TWO_DIMENSIONAL_SIDES;     )
_3D(auto side_list = THREE_DIMENSIONAL_SIDES;   )

    for(auto cell: *m_cells) {
        for(auto side: side_list) {
            auto face = cell->faces(side)->single();

            if (face->neighbor_exist(cell)) {
                // Сосед либо с этого процесса, либо удаленный
                auto neighbor = face->neighbor(cell);
                if (neighbor->is_local()) {
                    continue;
                }

                auto nei_rank = neighbor->owner();

                // Добавим ячейку в граничный список
                if (borders.count(nei_rank) < 1) {
                    borders[nei_rank] = NodeList::create(nei_rank, ListRole::BORDER);
                    ghosts[nei_rank] = NodeList::create(nei_rank, ListRole::GHOST);
                }
                borders[nei_rank]->push_back(cell);
                ghosts[nei_rank]->push_back(neighbor);
            }
            else {
                // Сосед отсутствует - фейковая ячейка
                auto neighbor = create_neighbor(cell, side);

                connect_cells(side, cell, neighbor);
                int nei_rank = m_decomposition->rank(neighbor->center());

                if (nei_rank < 0) {
                    // Установим граничное условие
                    FaceFlag flag = face_flags[to_int(side)];

                    cell->faces(side)->single()->set_flag(flag);

                    if (flag == FaceFlag::ZOE) {
                        neighbor->set_data(cell->data_holder());
                    }
                    if (flag == FaceFlag::PERIODIC) {
                        DataHolder::Ptr real_data;
                        static size_t Nx = ((RectangleGenerator*)m_geometry.get())->global_nx() - 1;
                        if (face->normal(neighbor)[0] > 0.0) {
                            real_data = m_base_cells[cell->z_base() - Nx]->data_holder();
                        } else {
                            real_data = m_base_cells[cell->z_base() + Nx]->data_holder();
                        }
                        neighbor->set_data(real_data);
                    }

                    m_fake->push_back(neighbor);
                }
                else {
                    throw runtime_error("Non fake neighbor doest not exists");
                }
            }
        }
    }

    for(auto& pair: borders) {
        m_border.emplace_back(BorderControl::create(pair.second));
    }
    for(auto& pair: ghosts) {
        m_ghosts.emplace_back(GhostControl::create(pair.second));
    }
}

void Mesh::connect_cells(Side side, Cell::Ref cell, Cell::Ref neighbor) {
    auto face = cell->faces(side)->single();

    if (face->outward_normal().dot(neighbor->center() - cell->center()) > 0.0) {
        cell->faces(side)->single()->set_inward_cell(cell);
        cell->faces(side)->single()->set_outward_cell(neighbor);
    }
    else {
        cell->faces(side)->single()->set_outward_cell(cell);
        cell->faces(side)->single()->set_inward_cell(neighbor);
    }
}

Cell::Ptr Mesh::create_neighbor(Cell::Ref cell, Side side) {
    // formatter:off
    if (cell->faces(side)->is_complex()) {
        throw runtime_error("Mesh::create_neighbor(...) error: Attempt to create neighbor for complex face");
    }

    // Две грани, которые образуют новую ячейку
    auto face_1 = cell->faces(side)->single();
    auto face_2 = m_geometry->next_face_by_side(cell, side);

    static_vector<Face::Ptr, 2 * DIM> faces;
    Cell::Ptr neighbor;

    faces.resize(2 * DIM);

#if DIM2 /// TWO-DIMENSIONAL VERSION

    faces[to_int(opposite_side(side))] = face_1;
    faces[to_int(side)] = face_2;

    if (side == Side::LEFT) {

        faces[to_int(Side::BOTTOM)] = Face::create(face_1->vertex(0), face_2->vertex(0));
        faces[to_int(Side::TOP)]    = Face::create(face_1->vertex(1), face_2->vertex(1));

    } else if (side == Side::BOTTOM) {

        faces[to_int(Side::LEFT)]   = Face::create(face_2->vertex(1), face_1->vertex(1));
        faces[to_int(Side::RIGHT)]  = Face::create(face_2->vertex(0), face_1->vertex(0));

    } else if (side == Side::RIGHT) {

        faces[to_int(Side::BOTTOM)] = Face::create(face_2->vertex(0), face_1->vertex(0));
        faces[to_int(Side::TOP)]    = Face::create(face_2->vertex(1), face_1->vertex(1));

    } else if (side == Side::TOP) {

        faces[to_int(Side::LEFT)]   = Face::create(face_1->vertex(1), face_2->vertex(1));
        faces[to_int(Side::RIGHT)]  = Face::create(face_1->vertex(0), face_2->vertex(0));
    }

    neighbor = Cell::create(
            faces[to_int(Side::LEFT)],
            faces[to_int(Side::BOTTOM)],
            faces[to_int(Side::RIGHT)],
            faces[to_int(Side::TOP)]);

#else   /// THREE-DIMENSIONAL VERSION
    faces[to_int(opposite_side(side))] = face_1;
    faces[to_int(side)] = face_2;

    if (side == Side::LEFT) {

        faces[to_int(Side::BOTTOM)] = Face::create(face_2->vertex(0), face_2->vertex(3),
                                                   face_1->vertex(3), face_1->vertex(0));
        faces[to_int(Side::TOP)]    = Face::create(face_2->vertex(1), face_2->vertex(2),
                                                   face_1->vertex(2), face_1->vertex(1));
        faces[to_int(Side::BACK)]   = Face::create(face_2->vertex(0), face_1->vertex(0),
                                                   face_1->vertex(1), face_2->vertex(1));
        faces[to_int(Side::FRONT)]  = Face::create(face_2->vertex(3), face_1->vertex(3),
                                                   face_1->vertex(2), face_2->vertex(2));

    } else if (side == Side::RIGHT) {

        faces[to_int(Side::BOTTOM)] = Face::create(face_1->vertex(0), face_1->vertex(3),
                                                   face_2->vertex(3), face_2->vertex(0));
        faces[to_int(Side::TOP)]    = Face::create(face_1->vertex(1), face_1->vertex(2),
                                                   face_2->vertex(2), face_2->vertex(1));
        faces[to_int(Side::BACK)]   = Face::create(face_1->vertex(0), face_2->vertex(0),
                                                   face_2->vertex(1), face_1->vertex(1));
        faces[to_int(Side::FRONT)]  = Face::create(face_1->vertex(3), face_2->vertex(3),
                                                   face_2->vertex(2), face_1->vertex(2));

    } else if (side == Side::BOTTOM) {

        faces[to_int(Side::LEFT)]   = Face::create(face_2->vertex(0), face_1->vertex(0),
                                                   face_1->vertex(1), face_2->vertex(1));
        faces[to_int(Side::RIGHT)]  = Face::create(face_2->vertex(3), face_1->vertex(3),
                                                   face_1->vertex(2), face_2->vertex(2));
        faces[to_int(Side::BACK)]   = Face::create(face_2->vertex(0), face_2->vertex(3),
                                                   face_1->vertex(3), face_1->vertex(0));
        faces[to_int(Side::FRONT)]  = Face::create(face_2->vertex(1), face_2->vertex(2),
                                                   face_1->vertex(2), face_1->vertex(1));

    } else if (side == Side::TOP) {

        faces[to_int(Side::LEFT)]   = Face::create(face_1->vertex(0), face_2->vertex(0),
                                                   face_2->vertex(1), face_1->vertex(1));
        faces[to_int(Side::RIGHT)]  = Face::create(face_1->vertex(3), face_2->vertex(3),
                                                   face_2->vertex(2), face_1->vertex(2));
        faces[to_int(Side::BACK)]   = Face::create(face_1->vertex(0), face_1->vertex(3),
                                                   face_2->vertex(3), face_2->vertex(0));
        faces[to_int(Side::FRONT)]  = Face::create(face_1->vertex(1), face_1->vertex(2),
                                                   face_2->vertex(2), face_2->vertex(1));

    } else if (side == Side::BACK) {

        faces[to_int(Side::LEFT)]   = Face::create(face_2->vertex(0), face_2->vertex(3),
                                                   face_1->vertex(3), face_1->vertex(0));
        faces[to_int(Side::RIGHT)]  = Face::create(face_2->vertex(1), face_2->vertex(2),
                                                   face_1->vertex(2), face_1->vertex(1));
        faces[to_int(Side::BOTTOM)] = Face::create(face_2->vertex(0), face_1->vertex(0),
                                                   face_1->vertex(1), face_2->vertex(1));
        faces[to_int(Side::TOP)]    = Face::create(face_2->vertex(3), face_1->vertex(3),
                                                   face_1->vertex(2), face_2->vertex(2));

    } else if (side == Side::FRONT) {

        faces[to_int(Side::LEFT)]   = Face::create(face_1->vertex(0), face_1->vertex(3),
                                                   face_2->vertex(3), face_2->vertex(0));
        faces[to_int(Side::RIGHT)]  = Face::create(face_1->vertex(1), face_1->vertex(2),
                                                   face_2->vertex(2), face_2->vertex(1));
        faces[to_int(Side::BOTTOM)] = Face::create(face_1->vertex(0), face_2->vertex(0),
                                                   face_2->vertex(1), face_1->vertex(1));
        faces[to_int(Side::TOP)]    = Face::create(face_1->vertex(3), face_2->vertex(3),
                                                   face_2->vertex(2), face_1->vertex(2));

    }

    neighbor = Cell::create(
            faces[to_int(Side::LEFT)],
            faces[to_int(Side::BOTTOM)],
            faces[to_int(Side::RIGHT)],
            faces[to_int(Side::TOP)],
            faces[to_int(Side::FRONT)],
            faces[to_int(Side::BACK)]
    );
#endif

    neighbor->data_holder()->resize(cell->data_holder()->size());

    for (auto &face: faces) {
        if (face->outward_normal().dot(face->center() - neighbor->center()) > 0.0) {
            face->set_inward_cell(neighbor);
        } else {
            face->set_outward_cell(neighbor);
        }
    }
    return neighbor;
}

void Mesh::rebuild_ghosts() {
    for (auto &border: m_border) {
        border->rebuild();
    }
    for (auto &ghost:  m_ghosts) {
        ghost->rebuild();
    }
}

void Mesh::remove_ghost_connections() {
    throw std::runtime_error("Invalid 3D code");
    /*
    // Список базовых ячеек для гост-слоев
    NodeList::Ptr ghosts_base = NodeList::create(mpi::rank());
    for (auto &ghost: m_ghosts) {
        for (auto cell: *ghost->cells()) {
            ghosts_base->push_back(Cell::origin(cell));
        }
    }

    // Гост-ячейки с поверхности (могут быть скрытыми и не содержаться в m_ghosts)
    NodeList::Ptr ghosts_surf = NodeList::create(mpi::rank());
    for (auto cell: *ghosts_base) {
        vector<Cell::Ptr> leaf;
        Cell::leaf_children(cell, leaf);

        for (auto &lch: leaf) {
            ghosts_surf->push_back(lch);
        }
    }

    for (auto cell: *ghosts_surf) {
        for (auto side: TWO_DIMENSIONAL_SIDES) {
            auto faces = cell->faces(side);
            if (faces->is_simple()) {
                auto face = faces->single();
                if (face->neighbor_exist(cell)) {
                    auto nei = face->neighbor(cell);

                    if (nei->is_remote()) {
                        face->set_inward_cell(nullptr);
                        face->set_outward_cell(nullptr);
                    }
                }
                else {
                    face->set_inward_cell(nullptr);
                    face->set_outward_cell(nullptr);
                }
            }
            if (faces->is_complex()) {
                auto face_1 = faces->at(0);
                auto face_2 = faces->at(1);

                auto ex_1 = face_1->neighbor_exist(cell);
                auto ex_2 = face_2->neighbor_exist(cell);

                // Два соседа
                if (!ex_1 && !ex_2) {
                    auto new_face = Face::create(face_1->vertex(0), face_2->vertex(1));
                    cell->set_faces(side, FacesList::create(new_face));
                }
            }
        }
    }
     */
}

void Mesh::shuffle_and_split() {
    m_cells->shuffle();
    m_cells->split(m_n_chunks);
}

size_t Mesh::z_curve_offset_2d(size_t level) const {
    static vector<size_t> cells;
    if (cells.empty()) {
        cells = vector<size_t>(m_max_level + 1);

        cells[0] = m_geometry->n_cells();
        for(uint l = 1; l < m_max_level + 1; ++l) {
            cells[l] = cells[l - 1] + math::pow4(l) * cells[0];
        }
    }
    return cells[level];
}

size_t Mesh::z_curve_offset_3d(size_t level) const {
    static vector<size_t> cells;
    if (cells.empty()) {
        cells = vector<size_t>(m_max_level + 1);

        cells[0] = m_geometry->n_cells();
        for(uint l = 1; l < m_max_level + 1; ++l) {
            cells[l] = cells[l - 1] + math::pow2(3 * l) * cells[0];
        }
    }
    return cells[level];
}

// ============================================================================



// ============================================================================
//                              МЕТОДЫ ДОСТУПА
// ============================================================================

size_t Mesh::max_level() {
    return m_max_level;
}

NodeList::Ref Mesh::cells() {
    return m_cells;
}

uint Mesh::n_chunks() const {
    return m_cells->n_parts();
}

NodeList::Part Mesh::all_chunks() {
    return m_cells->all_parts();
}

const GeometryGenerator* Mesh::geometry() const {
    return m_geometry.get();
}

double Mesh::map(const std::function<void(const NodeList::Part &)>& func) {
    auto begin = std::chrono::steady_clock::now();
    m_n_chunks = m_cells->n_parts();
    if (m_n_chunks <= 1) {
        func(m_cells->all_parts());
    }
    else {
        vector<thread> workers(m_n_chunks);
        for (uint i = 0; i < m_n_chunks; ++i) {
            workers[i] = thread(func, m_cells->part(i));
        }
        for (auto &worker: workers)
            worker.join();
    }
    auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(end - begin).count();
}

void Mesh::extend_cell_data(size_t cell_data_size) {
    for(auto cell: *m_cells) {
        cell->data_holder()->resize(cell_data_size);
    }
    for(auto cell: *m_fake) {
        cell->data_holder()->resize(cell_data_size);
    }
    for(auto& ghost: m_ghosts) {
        for(auto cell: *ghost->cells()) {
            cell->data_holder()->resize(cell_data_size);
        }
    }
}

// ============================================================================



// ============================================================================
//                                ОТЛАДКА
// ============================================================================

void Mesh::dbg_write_ghost_connectivity(size_t step) {
    for(size_t gid = 0; gid < m_ghosts.size(); ++gid) {
        auto path = "ghost_connectivity_" + to_string(mpi::rank()) + "_" + to_string(gid) + "_" ;

        size_t write_number = 0;

        vector<Cell::Ptr> to_write;
        for(auto cell: *m_ghosts[gid]->cells()) {
            to_write.push_back(cell);
            cell->set_id(0);
        }

        for(auto& cell: to_write) {
            auto faces = cell->faces_list();

            cell->set_id(faces.size());

            for(auto face: faces) {
                if(face->neighbor_exist(cell)) {
                    face->neighbor(cell)->set_id(1);
                }
                else {
                    cell->set_id(42);
                }
            }

            VtkWriter::write(*m_cells, "mesh_" + path + to_string(write_number) + ".vtu");
            VtkWriter::write(*m_ghosts[gid]->cells(), path + to_string(write_number) + ".vtu");

            write_number++;

            for(auto face: faces) {
                if(face->neighbor_exist(cell)) {
                    face->neighbor(cell)->set_id(0);
                }
            }

            cell->set_id(0);
        }
    }
}

void Mesh::dbg_write_connectivity(size_t step) {
    auto path = "connectivity_" + to_string(mpi::rank()) + "_" + to_string(step) + "_" ;

    size_t write_number = 0;

    vector<Cell::Ptr> to_write;

    for(auto node: *m_cells) {
        auto cell = std::static_pointer_cast<Cell>(node);
        to_write.push_back(cell);
        cell->set_adaptation_flag(static_cast<AdaptationFlag>(0));
    }

    static int cell_id = 0;

    for(auto& cell: to_write) {
        auto faces = cell->faces_list();

        cell->set_adaptation_flag(static_cast<AdaptationFlag>(faces.size()));


        auto neighbors_exist = cell->neighbors().size();
        if(neighbors_exist == faces.size()) {
            cell->set_adaptation_flag(static_cast<AdaptationFlag>(neighbors_exist));
        }
        else {
            cell->set_adaptation_flag(static_cast<AdaptationFlag>(20 + neighbors_exist));

            auto nl = NodeList::create(cell->owner());
            nl->push_back(cell);

            for(auto& n: cell->neighbors()) {
                nl->push_back(n);
            }

            VtkWriter::write(*nl, "nn_" + to_string(mpi::rank()) + "_" + to_string(cell_id++) + ".vtu");

            nl->clear();
        }

        for(auto& face: faces) {
            if (face->neighbor_exist(cell)) {
                face->neighbor(cell)->set_adaptation_flag(static_cast<AdaptationFlag>(1));
            }
        }

        VtkWriter::write(*m_cells, path + to_string(++write_number) + ".vtu");

        for(auto& face: faces) {
            if(face->neighbor_exist(cell)) {
                reset_adaptation_flag(face->neighbor(cell));
            }
        }

        reset_adaptation_flag(cell);
    }
}

void Mesh::dbg_write(Problem* problem, const string& prefix) {
    static int step = 0;

    auto to_str = [](int x) -> string {
        string num = to_string(x);
        if (num.length() < 3) {
            return string(3 - num.length(), '0') + num;
        }
        else {
            return num;
        }
    };

    string dir = "debug/";
    string file = "dbg_" + to_str(step);
    if (!prefix.empty()) {
        file += "_" + prefix;
    }


    if (mpi::is_master()) {
        std::ofstream pvd(dir + file + ".pvd");

        pvd << "<?xml version=\"1.0\"?>\n";
        pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        pvd << "    <Collection>\n";

        for (int r = 0; r < mpi::size(); ++r) {
            // main
            pvd << "        <DataSet timestep=\"0.0\""
                << " part=\"" << r << "\""
                << " file=\"" << file << "_A.pt" << r << ".vtu\"/>\n";

            // fake
            pvd << "        <DataSet timestep=\"0.0\""
                << " part=\"" << mpi::size() + r << "\""
                    << " file=\"" << file << "_B(fk).pt" << r << ".vtu\"/>\n";
        }

        pvd << "    </Collection>\n";
        pvd << "</VTKFile>\n";

        pvd.close();
    }

    string rank = mpi::srank();

    auto filename_main = dir + file + "_A.pt" + rank + ".vtu";
    auto filename_fake = dir + file + "_B(fk).pt" + rank + ".vtu";
    auto filename_border = dir + file + "_C(bd).pt" + rank + ".";
    auto filename_ghost = dir + file + "_D(gh).pt" + rank + ".";
    auto filename_sghost = dir + file + "_E(sgh).pt" + rank + ".vtu";
    auto filename_base = dir + file + "_F.pt" + rank + ".vtu";

    VtkWriter::write(*m_cells, filename_main, problem);
    if (m_fake != nullptr) {
        VtkWriter::write(*m_fake, filename_fake, nullptr);
    }

    for (auto &border: m_border) {
        VtkWriter::write(*border->cells(), filename_border + to_string(border->rank()) + ".vtu", problem);
    }
    for (auto &ghost:  m_ghosts) {
        VtkWriter::write(*ghost->cells(), filename_ghost + to_string(ghost->rank()) + ".vtu", problem);
    }
        
    // Записываем супер-гост ячейки
    NodeList::Ptr ghosts_base = NodeList::create(mpi::rank());
    for(auto& ghost: m_ghosts) {
        for(auto cell: *ghost->cells()) {
            ghosts_base->push_back(Cell::origin(cell));
        }
    }
    NodeList ghosts_surf(mpi::rank());
    for(auto cell: *ghosts_base) {
        vector<Cell::Ptr> leaf;
        Cell::leaf_children(cell, leaf);

        for(auto &lch: leaf) {
            ghosts_surf.push_back(lch);
        }
    }
    VtkWriter::write(ghosts_surf, filename_sghost);

    // Базовые ячейки
    NodeList base_cells(mpi::rank());
    for(auto cell: m_base_cells) {
        base_cells.push_back(cell.second);
    }
    VtkWriter::write(base_cells, filename_base);

    step++;

}
void Mesh::dbg_write_seq() {
    static int step = 0;
    int gc = 0;
    int bc = 0;

    string file = "dbg_seq_" + to_string(mpi::rank()) + "_";
    auto filename_ghost  = file + "ghost_";
    auto filename_border = file + "border_";

    for(auto& ghost: m_ghosts) {
        vector<Cell::Ptr> cells;

        for(auto node: *ghost->cells()) {
            auto cell = std::static_pointer_cast<Cell>(node);
            cell->set_id(0); cells.push_back(cell);
        }

        int i = 0;
        for(auto& cell: cells) {
            cell->set_id(1);

            VtkWriter::write(*ghost->cells(), filename_ghost + to_string(gc) + "_" + to_string(step) + "_" + to_string(i++) + ".vtu");

            cell->set_id(0);
        }

        gc++;
    }

    for(auto& border: m_border) {
        vector<Cell::Ptr> cells;

        for(auto node: *border->cells()) {
            auto cell = std::static_pointer_cast<Cell>(node);
            cell->set_id(0); cells.push_back(cell);
        }

        int i = 0;
        for(auto& cell: cells) {
            cell->set_id(1);

            VtkWriter::write(*border->cells(), filename_border + to_string(bc) + "_" + to_string(step) + "_" + to_string(i++) + ".vtu");

            cell->set_id(0);
        }

        bc++;
    }

    // writer->write(m_cells, file + "base");

    // for(auto& border: m_border) {
    //     writer->write(border->cells(), filename_border + to_string(bc) + "_" + to_string(solution_step) + ".vtu");
    //     bc++;
    // }

    step++;
}

void Mesh::dbg_check_duplicates() {
    vector<vector<size_t>> coords;

    for(auto node: *m_cells) {
        auto cell = std::static_pointer_cast<Cell>(node);

        for(auto& coord: coords) {
            if(coord[0] == cell->level() &&
               coord[1] == cell->z()) {
                throw runtime_error("Error! Cell duplicated in m_cells!");
            }
        }

        coords.push_back({cell->level(), cell->z()});
    }
}

void Mesh::print_info(const std::string& tab) {
    std::hash<double> d_hash;

    std::cout << tab << "Sector info:\n";
    m_geometry->print_info(tab + "  ");

    std::cout << tab << "Mesh info:\n";
    std::cout << tab << "  n: " << m_geometry->n_cells() << "\n";

    size_t z = 25;

    std::cout << tab << "Cell (" << z << "):\n";
    auto cell = base_cell(z);
    auto center = cell->center_2d();
    double r = cell->center_radius();
    double alpha = cell->center_angle();

    std::cout << tab << "  center_2d (cartesian): "
              << center[0] << ", "
              << center[1] << "\n";

    std::cout << tab << "  center_2d (cartesian): "
              << d_hash(center[0]) << ", "
              << d_hash(center[1]) << "\n";

    std::cout << tab << "  center_2d (polar):     "
              << d_hash(r) << ", "
              << d_hash(alpha) << "\n";

    std::cout << tab << "  cell data hash: " << cell->data_holder()->hash() << "\n";

    std::cout << tab << "full mesh hash: " << hash() << "\n";
}

size_t Mesh::hash() {
    size_t full_hash = 0;
    for(auto cell: *m_cells) {
        full_hash += cell->id() * cell->data_holder()->hash();
    }
    return full_hash;
}

// ============================================================================