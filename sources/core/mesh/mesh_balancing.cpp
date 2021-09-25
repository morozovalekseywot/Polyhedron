//
// Created by 159-mrv on 6/22/18.
//

#include <control/configuration.h>
#ifdef ENABLE_MPI
#include <control/mpi_commutator.h>
#endif

#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/data_holder.h>
#include <core/layer/border_control.h>
#include <core/layer/ghost_control.h>
#include <core/mesh/mesh.h>
#include <core/generator/structured_decomposition.h>
#include <core/generator/sector_generator.h>
#include <core/generator/decomposition.h>

#include <io/vtk/vtk_writer.h>
#include <problems/problem.h>

#include <core/cell/side.h>
#include <utils/math.h>
#include <utils/partitioning/voronoi_diagram.h>
#include <utils/memory/kd_tree.h>
#include <utils/partitioning/base_splitter_1d.h>
#include <utils/partitioning/base_splitter_2d.h>
#include <utils/partitioning/rcb_splitter.h>
#include <core/mesh/partitioner.h>
#include <utils/geom.h>
#include <core/generator/orth_decomposition.h>

vector<std::set<Cell::Ptr, CompareByIndex>> Mesh::compute_cells_to_send() const {
    vector<std::set<Cell::Ptr, CompareByIndex>> cells_to_send(mpi::usize());

    // ВНИМАНИЕ! Отправлять будем только собственные ячейки,
    // то есть только ячейки из m_fake и m_cells.
    // Отправлять гостовые ячейки нет необходимости.
    for (auto cell: *m_cells) {
        // Ячейка чужого процесса, надо отправить
        if (cell->is_remote()) {
            cells_to_send[cell->owner()].insert(cell);
        }

        // Данная ячейка граничит с другим процессом,
        // надо её отправить, может оказаться ghost ячейкой
        for (const auto &nei: cell->neighbors()) {
            if (nei->is_remote()) {
                cells_to_send[nei->owner()].insert(cell);
            }
        }
    }
    for (auto cell: *m_fake) {
        // Ячейка чужого процесса, надо отправить
        if (cell->is_remote()) {
            cells_to_send[cell->owner()].insert(cell);
        }

        // Данная ячейка граничит с другим процессом,
        // надо её отправить, может оказаться ghost ячейкой
        for (const auto &nei: cell->neighbors()) {
            if (nei->is_remote()) {
                cells_to_send[nei->owner()].insert(cell);
            }
        }
    }

    cells_to_send[mpi::rank()].clear();
    return cells_to_send;
}

#ifdef ENABLE_MPI
vector<ClearVertex> get_new_vertices(const vector<std::set<Cell::Ptr, CompareByIndex>>& cells_to_send, const MpiCommutator& comm) {
    vector<vector<ClearVertex>> clear_vertices_to_send(mpi::usize());
    for (int r = 0; r < mpi::size(); ++r) {
        std::set<Vertex::Ptr, CompareByIndex> vertices_to_send;

        for (const auto &cell: cells_to_send[r]) {
            for (auto& face: cell->faces_list()) {
                for (auto& v: face->vertices()) {
                    vertices_to_send.insert(v);
                }
            }
        }

        for(const auto& vertex: vertices_to_send) {
            clear_vertices_to_send[r].emplace_back(ClearVertex(vertex));
        }
    }

    return comm.send_to_common_buffer(clear_vertices_to_send);
}

vector<ClearFace> get_new_faces(const vector<std::set<Cell::Ptr, CompareByIndex>>& cells_to_send, const MpiCommutator& comm) {
    vector<ClearVertex> new_faces;

    // Запакуем в массив
    vector<vector<ClearFace>> clear_faces_to_send(mpi::usize());

    for(int r = 0; r < mpi::size(); ++r) {
        std::set<Face::Ptr, CompareByIndex> faces_to_send;

        for (const auto &cell: cells_to_send[r]) {
            for (const auto &face: cell->faces_list()) {
                faces_to_send.insert(face);
            }
        }

        for(const auto& face: faces_to_send) {
            clear_faces_to_send[r].emplace_back(ClearFace(face));
        }
    }

    return comm.send_to_common_buffer(clear_faces_to_send);
}

vector<ClearCell> get_new_cells(const NodeList::Ptr& fake_list,
                                const vector<std::set<Cell::Ptr, CompareByIndex>>& cells_to_send,
                                const MpiCommutator& comm) {
    vector<ClearCell> new_cells;

    vector<vector<ClearCell>> clear_cells_to_send(mpi::usize());

    for(int r = 0; r < mpi::size(); ++r) {
        clear_cells_to_send[r].resize(cells_to_send[r].size());

        size_t counter = 0;
        for(const auto& cell: cells_to_send[r]) {
            clear_cells_to_send[r][counter++] = ClearCell(cell, fake_list->contain(cell));
        }
    }

    return comm.send_to_common_buffer(clear_cells_to_send);
}

vector<double> get_new_data(const vector<std::set<Cell::Ptr, CompareByIndex>>& cells_to_send,
                            const MpiCommutator& comm) {
    vector<vector<double>> cell_data_to_send(mpi::usize());

    size_t cell_data_size = 0;
    for(int r = 0; r < mpi::size(); ++r) {
        if (!cells_to_send[r].empty()) {
            auto cell = *cells_to_send[r].begin();
            cell_data_size = cell->data_holder()->size();
        }
    }

    for(int r = 0; r < mpi::size(); ++r) {
        cell_data_to_send[r].resize(cell_data_size * cells_to_send[r].size());

        size_t counter = 0;
        for(const auto& cell: cells_to_send[r]) {
            cell->data_holder()->serialize(&cell_data_to_send[r][counter]);
            counter += cell_data_size;
        }
    }

    return comm.send_to_common_buffer(cell_data_to_send);
}
#endif

void Mesh::load_balance(double proc_workload) {
    if (!m_balancing_options.use_balancing()) {
        return;
    }

    if (proc_workload == 0.0) {
        proc_workload = m_cells->size();
    }

    double max_workload = mpi::max(proc_workload);
    double min_workload = mpi::min(proc_workload);
    double imbalance = (max_workload - min_workload) / (max_workload + min_workload);

    if (imbalance < 0.01) {
        return;
    }

    //dbg_write(nullptr, "lb_bef");

    auto prev_n_chunks = m_cells->n_parts();

    // Глобальная индексация примитивов, требуется в обменных операциях
    {
        std::map<IndexedObject::index_type, Vertex::Ptr> temp_verts;
        std::map<IndexedObject::index_type, Face::Ptr> temp_faces;
        std::map<IndexedObject::index_type, Cell::Ptr> temp_cells;

        global_indexing(temp_verts, temp_faces, temp_cells);
    }

    // Определяем следующего владельца ячеек
    compute_next_owners(proc_workload);

    //dbg_write(nullptr, "comp");

    // Отправляем и получаем новые примитивы в запакованном виде
    vector<ClearVertex> new_vertices;
    vector<ClearFace> new_faces;
    vector<ClearCell> new_cells;
    vector<double> new_data;
    {
        // Ячейки на отправку
        auto cells_to_send = compute_cells_to_send();

        // Определим потоки ячеек
        vector<bool> need_to_send(mpi::usize(), false);
        for (int r = 0; r < mpi::size(); ++r) {
            need_to_send[r] = !cells_to_send[r].empty();
        }
#ifdef ENABLE_MPI
        MpiCommutator comm(need_to_send);
        new_vertices = get_new_vertices(cells_to_send, comm);
        new_faces = get_new_faces(cells_to_send, comm);
        new_cells = get_new_cells(m_fake, cells_to_send, comm);
        new_data = get_new_data(cells_to_send, comm);
#endif
    }


    // ========================================================================
    // Теперь необходимо всё вычистить
    // ========================================================================

    // Соберем нужные ячейки
    std::map<IndexedObject::index_type, Vertex::Ptr> old_verts;
    std::map<IndexedObject::index_type, Face::Ptr> old_faces;
    std::map<IndexedObject::index_type, Cell::Ptr> old_cells;

    for(auto cell: *m_cells) {
        if (cell->is_remote()) {
            bool useless = true;
            for(auto& nei: cell->neighbors()) {
                if (nei->is_local()) {
                    useless = false;
                    break;
                }
            }
            if (useless) {
                continue;
            }
        }

        old_cells[cell->get_index()] = cell;
        for(auto& face: cell->faces_list()) {
            old_faces[face->get_index()] = face;
            for(auto& v: face->vertices()) {
                old_verts[v->get_index()] = v;
            }
        }
    }
    for(auto cell: *m_fake) {
        if (cell->is_remote()) {
            bool useless = true;
            for(auto& nei: cell->neighbors()) {
                if (nei->is_local()) {
                    useless = false;
                    break;
                }
            }
            if (useless) {
                continue;
            }
        }

        old_cells[cell->get_index()] = cell;
        for(auto& face: cell->faces_list()) {
            old_faces[face->get_index()] = face;
            for(auto& v: face->vertices()) {
                old_verts[v->get_index()] = v;
            }
        }
    }

    // Удаляем ghost слой
    for(auto & ghost: m_ghosts) {
        for(auto cell: *ghost->cells()) {
            m_base_cells.erase(cell->z_base());
        }
        ghost->cells()->clear();
    }
    m_ghosts.clear();

    // Удаляем border слой
    for(auto& border: m_border) {
        border->cells()->clear();
    }
    m_border.clear();

    // Чистим m_cells
    for(auto it = m_cells->begin(); it != m_cells->end(); ) {
        auto cell = *it;
        ++it;
        if (cell->is_remote()) {
            cell->remove_from_all_lists();
        }
    }

    // Чистим m_fake
    for(auto it = m_fake->begin(); it != m_fake->end(); ) {
        auto cell = *it;
        ++it;
        if (cell->is_remote()) {
            cell->remove_from_all_lists();
        }
    }


    // ========================================================================
    // Распакуем чистые примитивы в нормальные
    // ========================================================================

    size_t cell_data_size = m_cells->begin()->data_holder()->size();
    vector<Cell::Ptr> for_hierarchy_restore;

    // parse all_verts
    for (const auto &v: new_vertices) {
        auto it = old_verts.find(v.index);
        // получили новую вершину
        if (it == old_verts.end()) {
            auto vertex = Vertex::create(v.x, v.y);
            vertex->set_index(v.index);

            old_verts[v.index] = vertex;
        }
    }
    new_vertices.clear();

    // parse all_faces
    for (const auto &f: new_faces) {
        auto fit = old_faces.find(f.index);

        // сторона существует, проехали
        if (fit != old_faces.end()) {
            continue;
        }

        static_vector<Vertex::Ptr, 4> vit;

        for(int i = 0; i < 4; ++i) {
            if (f.vertex[i] != IndexedObject::undefined_index) {
                auto vi = old_verts.find(f.vertex[i]);
                if (vi != old_verts.end()) {
                    vit.push_back(vi->second);
                }
            }
        }
        if (vit.size() == 2) {
            Face::Ptr face = Face::create(vit[0], vit[1]);
            face->set_index(f.index);
            old_faces[f.index] = face;
        }
        else {
            throw runtime_error("Unknown vertex on face");
        }
    }

    // parse all_cells
    for (size_t i = 0; i < new_cells.size(); ++i) {
        auto c = new_cells[i];

        // уже есть такая ячейка
        if (old_cells.find(c.index) != old_cells.end()) { continue; }

        array<FacesList::Ptr, 4> cell_faces;

        for (int s = 0; s < 4; ++s) {
            FacesList::Ptr side_faces;

            Face::Ptr face_1 = old_faces[c.faces[2 * s]];

            if (c.faces[2 * s + 1] != IndexedObject::undefined_index) {
                Face::Ptr face_2 = old_faces[c.faces[2 * s + 1]];
                side_faces = FacesList::create(face_1, face_2);
            }
            else {
                side_faces = FacesList::create(face_1);
            }

            cell_faces[s] = side_faces;
        }

        Cell::Ptr cell = Cell::create(cell_faces[to_int(Side::LEFT)],
                                      cell_faces[to_int(Side::BOTTOM)],
                                      cell_faces[to_int(Side::RIGHT)],
                                      cell_faces[to_int(Side::TOP)]);
        cell->set_id(c.id);
        cell->set_z(c.z);
        cell->set_level(c.level);
        cell->set_index(c.index);
        cell->set_owner(c.owner);

        cell->data_holder()->resize(cell_data_size);
        cell->data_holder()->deserialize(&new_data[cell_data_size * i]);

        old_cells[cell->get_index()] = cell;
        if (c.fake) {
            m_fake->push_back(cell);
        } else {
            if (cell->owner() == mpi::rank()) {
                m_cells->push_back(cell);
            }
            for_hierarchy_restore.push_back(cell);
        }
    }
    new_cells.clear();

    // снова парсим грани и восстанавливаем связи
    for (const auto &f: new_faces) {
        auto face = old_faces[f.index];

        if (f.inward_cell != IndexedObject::undefined_index) {
            auto it = old_cells.find(f.inward_cell);
            if (it != old_cells.end()) {
                face->set_inward_cell(it->second);
            }
        }

        if (f.outward_cell != IndexedObject::undefined_index) {
            auto it = old_cells.find(f.outward_cell);
            if (it != old_cells.end()) {
                face->set_outward_cell(it->second);
            }
        }

        face->set_flag(to_face_flag(f.flag));
    }
    new_faces.clear();

    restore_hierarchy(for_hierarchy_restore);

    for (auto &border: m_border) { border->cells()->clear(); }
    for (auto &ghost:  m_ghosts) { ghost->cells()->clear(); }

    m_border.clear();
    m_ghosts.clear();

    for (auto cell: *m_cells) {
        if (cell->owner() != mpi::rank()) {
            throw runtime_error("Error! Cell with non local owner!");
        }

        for (auto &neighbor: cell->neighbors()) {
            if (m_cells->contain(neighbor)) { continue; }
            if (m_fake->contain(neighbor)) { continue; }

            auto bid = add_neighbor_if_not_exist(neighbor->owner());

            if (neighbor->is_local()) {
                throw runtime_error("Error: Ghost cell with local owner");
            }

            push_cell_if_not_equal(cell, m_border[bid]->cells());
            push_cell_if_not_equal(neighbor, m_ghosts[bid]->cells());
        }
    }

    m_cells->split(prev_n_chunks);

    for (auto &border: m_border) { border->rebuild(); }
    for (auto &ghost:  m_ghosts) { ghost->rebuild(); }

    remove_ghost_connections();

    // На данный момент у гост-слоя может оставаться полная иехрархия
    // которая включает ячейки максимального уровня для всех сиблингов гост слоя
    clean_ghosts();

    //dbg_write(nullptr, "lb_aft");
}

void Mesh::compute_next_owners(double workload) {
    m_decomposition->set_workload(workload);
    m_decomposition->sync_data();
    m_decomposition->update();

    auto get_rank = [this](Cell::Ref cell) -> int {
        auto grandpa = Cell::origin(cell);
        return m_decomposition->rank(grandpa->center());
    };

    // Следующий владелец для ячеек сетки
    for (auto cell: *m_cells) {
        cell->set_owner(get_rank(cell));
    }

    // Следующий владелец для ghost ячеек
    for (auto &ghost: m_ghosts) {
        for (auto cell: *ghost->cells()) {
            cell->set_owner(get_rank(cell));
        }
    }

#if 0 // Voronoi
    vector<double> workloads = mpi::all_gather(workload);
    static bool first_time = true;
    static VoronoiDiagram vd(mpi::size());
    if (first_time) {
        vector<Vector2d> gens = vector<Vector2d>(mpi::usize());

        double x = 0.0;
        double y = 0.0;

        for(auto cell: *m_cells) {
            x += cell->center_2d()[0];
            y += cell->center_2d()[1];
        }

        x /= m_cells->size();
        y /= m_cells->size();

        srand(mpi::rank());

        double amp = 1.0e-9;
        //x += amp * (2.0 * rand() - 0.5);
        y += amp * (2.0 * rand() - 0.0);

        vector<double> xs = mpi::all_gather(x);
        vector<double> ys = mpi::all_gather(y);

        for (int i = 0; i < mpi::size(); ++i) {
            Vector2d gen = {xs[i], ys[i]};
            vd.set_generator(i, gen);
        }
        vd.update();

        first_time = false;
    }
    else {
        // Определяем центр масс и среднюю длину стороны ячейки
        vector<Vector2d> Rc = vector<Vector2d>(mpi::usize());
        double length;
        {
            Vector2d rc = Vector2d::Zero();
            double volume = 0.0;
            double rho_sum = 0.0;
            for (auto cell: *m_cells) {
                double rho = 1.0 / cell->volume();
                volume += cell->volume();
                rc += rho * cell->center_2d();
                rho_sum += rho;
            }
            volume = volume / m_cells->size();
            rc = rc / rho_sum;
            length = sqrt(volume);

            length = mpi::min(length);

            vector<double> rcx = mpi::all_gather(rc[0]);
            vector<double> rcy = mpi::all_gather(rc[1]);

            for(int i = 0; i < mpi::size(); ++i) {
                Rc[i] = {rcx[i], rcy[i]};
            }
        }

        // Определяем новое положение
        Vector2d gen = vd.get_generator(mpi::rank());
        Vector2d dg = Vector2d::Zero();
        Vector2d DG1 = Vector2d::Zero();
        Vector2d DG2 = Vector2d::Zero();


        for(size_t k = 0; k < m_ghosts.size(); ++k) {
            double Lij = math::inf();
            for(auto cell: *m_ghosts[k]->cells()) {
                Lij = min(Lij, sqrt(cell->volume()));
            }
            for(auto cell: *m_border[k]->cells()) {
                Lij = min(Lij, sqrt(cell->volume()));
            }

            int i = mpi::rank();
            int j = m_ghosts[k]->rank();

            double diff = (workloads[j] - workloads[i]) /
                          (workloads[j] + workloads[i]);

            //if (fabs(diff) < 1.0e-3) {
            //    diff = 0.0;
            //}

            double Dsh = 1.0 * Lij;
            std::cout << "\tdiff(" << i << ", " << j << "): " << diff << "\n";
            //std::cout << "\tL(" << i << ", " << j << "): " << Lij << "\n";

            DG1 = Dsh * diff * (vd.get_generator(j) - gen).normalized();
        }

        DG2 = Rc[mpi::rank()] - gen;

        double theta = 0.05;

        dg = (1 - theta) * DG1 + theta * DG2;

        // Сглаживаем векторное поле dg
        {
            vector<double> dxs = mpi::all_gather(dg[0]);
            vector<double> dys = mpi::all_gather(dg[1]);

            double dx = 0.0;
            double dy = 0.0;
            double sum_w = 0.0;

            int i = mpi::rank();
            for (auto &ghost: m_ghosts) {
                int j = ghost->rank();

                double w_i = 1.0 / (vd.get_generator(i) - vd.get_generator(j)).squaredNorm();

                dx += w_i * dxs[j];
                dy += w_i * dys[j];

                sum_w += w_i;
            }

            dx /= sum_w;
            dy /= sum_w;

            double gamma = 0.0;

            dg[0] = (1.0 - gamma) * dxs[i] + gamma * dx;
            dg[1] = (1.0 - gamma) * dys[i] + gamma * dy;
         }

        /// Новые значения
        double new_x = gen[0] + dg[0];
        double new_y = gen[1] + dg[1];

        // Направим генератор внутрь области
        {
            double Rout = 11050;
            if (new_y < 0.0) {
                new_y = 0.0;
            }
            double r2 = new_x * new_x + new_y * new_y;
            double r = sqrt(r2);
            if (r > Rout) {
                new_x = r * new_x / Rout;
                new_y = r * new_y / Rout;
            }
        }


        vector<double> new_gens_x = mpi::all_gather(new_x);
        vector<double> new_gens_y = mpi::all_gather(new_y);

        for(int i = 0; i < mpi::size(); ++i) {
            Vector2d v = {new_gens_x[i], new_gens_y[i]};
            vd.set_generator(i, v);
        }

        for(int i = 0; i < mpi::size(); ++i) {
            if (mpi::rank() == i) {
                std::cout << "\tg("<< i << "): "
                          << vd.get_generator(i)[0] << " "
                          << vd.get_generator(i)[1] << "\n";

                std::cout << "\tDG1(" << i << "): " << DG1[0] << " " << DG1[1] << "\n";
                std::cout << "\tDG2(" << i << "): " << DG2[0] << " " << DG2[1] << "\n";
                std::cout << "\tdg("<< i << "): " << dg[0] << " " << dg[1] << "\n";
            }
            mpi::barrier();
        }

        vd.update();
    }

    // Следующий владелец для ячеек сетки
    for (auto cell: *m_cells) {
        cell->set_owner(vd.cell_index(Cell::origin(cell)->center_2d()));
    }

    // Следующий владелец для ghost ячеек
    for (auto &ghost: m_ghosts) {
        for (auto cell: *ghost->cells()) {
            cell->set_owner(vd.cell_index(Cell::origin(cell)->center_2d()));
        }
    }
    */
#endif // Voronoi


    // Следующий владелец для фейковых ячеек
    for (auto cell: *m_fake) {
        for (const auto &nei: cell->neighbors()) {
            // нефейковый сосед,
            // эта фейковая ячейка нужна процессу этого соседа
            if (!m_fake->contain(nei)) {
                cell->set_owner(nei->owner());
            }
        }
    }
}

void Mesh::restore_hierarchy(CellVector &cells) {
    for(auto& cell: cells) {
        if (cell->level() != 0)
            continue;

        if (!base_cell_exist_for(cell)) {
            insert_base_cell(cell);
        }
    }

    for(size_t level = 1; level <= max_level(); ++level) {
        for(auto& cell: cells) {
            if(cell->level() != level) {
                continue;
            }

            // Если базовая ячейка не существует - создадим её
            if(!base_cell_exist_for(cell)) {
                auto batya = m_geometry->create_base_cell(cell->z_base());
                insert_base_cell(batya);
            }

            insert_child(base_cell(cell->z_base()), cell);
        }
    }
}

void Mesh::insert_child(Cell::Ptr grandpa, Cell::Ref child) {
    grandpa->set_owner(child->owner());

    // К данному моменту мы гарантировано имеем корневую ячейку и корневую соседнюю ячейку
    // (при этом они связаны соседством), почему бы не начать их связывать
    while (grandpa->level() + 1 < child->level()) {
        // К данному моменту neighbor_root однозначно существует
        // Измельчим его и выберем детей, в направлении целевой ячейки
        if (grandpa->is_childless()) {
            split_cell(grandpa, true);
        }

        grandpa = grandpa->child_directed_to(child);
    }

    // Последняя корневая соседская ячейка на уровень меньше той которую необходимо
    // интегрировать. Если у этой ячейки отсутствуют дети - создадим их
    if (grandpa->is_childless()) {
        split_cell(grandpa, true);
    }

    // К данному моменту мы имеем neighbor_root, в который
    // необходимо внедрить соседнию ячейку
    size_t loc_z = child->z() % 4;
    Cell::set_child_by_local_id(grandpa, loc_z, child);
}

size_t Mesh::add_neighbor_if_not_exist(int rank) {
    auto it = std::find_if
            (
                    m_ghosts.begin(),
                    m_ghosts.end(),
                    [&rank](GhostControl::URef ghost){ return ghost->rank() == rank; }
            );

    if(it == m_ghosts.end()) {
        m_ghosts.emplace_back(GhostControl::create(NodeList::create(rank)));
        m_border.emplace_back(BorderControl::create(NodeList::create(rank)));

        return m_ghosts.size() - 1;
    }
    else {
        return static_cast<size_t>(std::distance(m_ghosts.begin(), it));
    }
}

void Mesh::push_cell_if_not_equal(Cell::Ref cell, NodeList::Ptr cells) {
    auto this_rank = mpi::rank();

    if(cells->contain(cell)) { return; }

    static int step = 0;

    bool duplicate = false;
    for(auto gc: *cells) {

        if(cell->equal_to(gc)) {
            if(cell == gc) { continue; }

            std::cerr << "\n\n\n DUPLICATE IN ";

            if(cells->role() == ListRole::GHOST) {
                std::cerr << "GL (" << gc->z() << " - " << gc->level() << ")! \n\n\n";
            }
            else {
                std::cerr << "BL (" << gc->z() << " - " << gc->level() << ")! \n\n\n";
            }

            // tst
            auto gd = NodeList::create(this_rank); gd->push_back(gc);
            auto cd = NodeList::create(this_rank); cd->push_back(cell);

            VtkWriter::write(*gd, "GL_" + to_string(this_rank) + "_" + to_string(step) + ".vtu");
            VtkWriter::write(*cd, "CL_" + to_string(this_rank) + "_" + to_string(step) + ".vtu");

            gd->clear();
            cd->clear();

            dbg_write_connectivity(0);

            throw runtime_error("Error! Find duplicate!");
            step++;
            // tst
        }
    }

    if(!duplicate) {
        cells->push_back(cell);
    }
}

bool Mesh::base_cell_exist_for(Cell::Ref cell) {
    return m_base_cells.find(cell->z_base()) != m_base_cells.end();
}

void Mesh::insert_base_cell(Cell::Ref cell) {
    if(cell->level() != 0) {
        throw runtime_error("Error! Trying to add base cell that has level > root_level");
    }

    m_base_cells[cell->z()] = cell;
}

Cell::Ptr Mesh::base_cell(size_t z) {
    return m_base_cells[z];
}

void Mesh::clean_ghosts() {
    // Список базовых ячеек для гост-слоев
    NodeList::Ptr ghosts_base = NodeList::create(mpi::rank());
    for(auto& ghost: m_ghosts) {
        for(auto cell: *ghost->cells()) {
            ghosts_base->push_back(Cell::origin(cell));
        }
    }

    // Гост-ячейки с поверхности (могут быть скрытыми и не содержаться в m_ghosts)
    NodeList::Ptr ghosts_surf = NodeList::create(mpi::rank());
    for(auto cell: *ghosts_base) {
        vector<Cell::Ptr> leaf;
        Cell::leaf_children(cell, leaf);

        for(auto &lch: leaf) {
            ghosts_surf->push_back(lch);
        }
    }

    /*
    for(auto cell: *ghosts_surf) {
        if (!m_cells->contain(cell)) {
            set_need_to_coarse(cell);
        }
    }

    for(auto& ghost: m_ghosts) {
        for(auto cell: *ghost->cells()) {
            reset_adaptation_flag(cell);
        }
    }

     */
    static size_t counter = 0;

    //VTKWriter::write(*ghosts_base, "debug/gbase_" + to_string(counter) + ".pt" + mpi::srank() + ".vtu");
    //VTKWriter::write(*ghosts_surf, "debug/gsurf_" + to_string(counter) + ".pt" + mpi::srank() + ".vtu");

    ++counter;
}