//
// Created by 159-mrv on 4/23/18.
//

#include <control/configuration.h>
#include <control/mpi_wrapper.h>

#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/side.h>
#include <core/cell/cell.h>
#include <core/cell/data_holder.h>
#include <core/layer/border_control.h>
#include <core/layer/ghost_control.h>
#include <core/mesh/mesh.h>
#include <core/generator/sector_generator.h>
#include <core/generator/real_sector_generator.h>

#include <io/vtk/vtk_writer.h>
#include <problems/problem.h>
#include <utils/geom.h>
#include <utils/stopwatch.h>

void Mesh::geometric_adaptation(AdaptationCriterion criterion) {
    mpi::cout << "    Geometric adaptation\n";

    if (criterion == AdaptationCriterion::NONE) {
        return;
    }

    // Каскад поддерживается только для областей типа сектор
    if (criterion == AdaptationCriterion::CASCADE &&
        m_geometry->type() == GeometryType::RECTANGLE) {
        criterion = AdaptationCriterion::NONE;
    }

    size_t max_level = m_max_level;
    // Даже если адаптаця не нужна, делаем на первом этапе
    // Это необходимо для нормального разрешения начальных условий
    if (criterion == AdaptationCriterion::NONE) {
        max_level = m_max_level / 2;
        criterion = AdaptationCriterion::ALL;
    }

    for(size_t level = 0; level < max_level; ++level) {
        // выставляем флаги
        if (criterion == AdaptationCriterion::ALL) {
            map([this](const NodeList::Part& cells) {
                for(auto cell: cells) {
                    if (cell->level() < m_max_level) {
                        cell->set_adaptation_flag(AdaptationFlag::SPLIT);
                    }
                }
            });
        }
        else if (criterion == AdaptationCriterion::CASCADE) {
            map([this](const NodeList::Part &cells) {
                for (auto cell: cells) {
                    set_need_to_coarse(cell);
                }
            });

            map([this](const NodeList::Part& cells) {
                for(auto cell: cells) {
                    size_t wanted_level = level_by_radius(cell->center_radius());
                    if (wanted_level > cell->level()) {
                        set_need_to_split(cell);
                    }
                    else {
                        if (wanted_level == cell->level()) {
                            reset_adaptation_flag(cell);
                        }
                    }
                }
            });
        }
        else {
            throw runtime_error("Unknown type of geometric adaptation");
        }

        // Синхронизуем флаги
        sync_adaptation();

        // Применяем флаги
        apply_flags(nullptr);

        // Перестраиваем ghost и border
        rebuild_ghosts();

        // Переразбиваем чанки
        shuffle_and_split();
    }
}

void Mesh::adaptation(Problem* problem) {
    if (!m_solver_adaptation) {
        return;
    }

    reset_flags();
    set_flags(problem);
    sync_adaptation();
    //check_flags();

    apply_flags(problem);
    rebuild_ghosts();
    shuffle_and_split();

/*
    static size_t adaptation_step = 0;
    bool debug = true;
    ++adaptation_step;

    Stopwatch sw1;
    reset_flags();
    if (debug) { dbg_write(problem, "reset"); }
    set_flags(problem);
    sw1.stop();

    if (debug) { dbg_write(problem, "set"); }


    Stopwatch sw2;
    sync_adaptation();
    sw2.stop();

    if (debug) { dbg_write(problem, "sync"); }


    Stopwatch sw3;
    apply_flags(problem);
    sw3.stop();

    if (debug) { dbg_write(problem, "apply"); }


    Stopwatch sw4;
    rebuild_ghosts();
    shuffle_and_split();
    sw4.stop();

    for (int r = 0; r < mpi::size(); ++r) {
        mpi::barrier();
        if (r == mpi::rank()) {
            std::cout << "rank " << r << "\n";
            std::cout << "  set:     " << sw1.milliseconds() << " ms\n";
            std::cout << "  sync:    " << sw2.milliseconds() << " ms\n";
            std::cout << "  apply:   " << sw3.milliseconds() << " ms\n";
            std::cout << "  rebuild: " << sw4.milliseconds() << " ms\n";
        }
    }
*/
}

AdaptationFlag Mesh::adaptation_criterion(Problem* problem, Cell::Ref cell) {
    AdaptationFlag sflag = problem->adaptation_criterion(cell);

    AdaptationFlag flag = AdaptationFlag::NONE;
    size_t level;
    switch (m_adaptation_type) {
        case AdaptationCriterion::NONE:
            flag = sflag;
            break;

        case AdaptationCriterion::ALL:
            flag = AdaptationFlag::SPLIT;
            break;

        case AdaptationCriterion::CASCADE:
            level = level_by_radius(cell->center_radius());

            if (cell->level() < level) {
                flag = AdaptationFlag::SPLIT;
            }
            else {
                if (cell->level() == level) {
                    flag = sflag != AdaptationFlag::COARSE ? sflag : AdaptationFlag::NONE;
                }
                else {
                    // cell->level() > level
                    // можно делать то, что хочет решатель
                    flag = sflag;
                }
            }
            break;

//        default:
//            flag = AdaptationFlag::COARSE;
            break;
    }

    if (flag == AdaptationFlag::SPLIT && cell->level() >= m_max_level) {
        return AdaptationFlag::NONE;
    }
    // Нельзя огрублять "дедов"
    if (flag == AdaptationFlag::COARSE && cell->parent()) {
        bool grandpa = false;
        for(auto& child: cell->parent()->children()) {
            if (!child->children().empty()) {
                grandpa = true;
                break;
            }
        }

        flag = grandpa ? AdaptationFlag::NONE : AdaptationFlag::COARSE;
    }
    return flag;
}

void Mesh::reset_flags() {
    map([this](const NodeList::Part &cells) {
        for (auto cell: cells) {
            set_need_to_coarse(cell);
        }
    });

    for(auto& ghost: m_ghosts) {
        for(auto cell: *ghost->cells()) {
            set_need_to_coarse(cell);
        }
    }
}

void Mesh::set_flags(Problem* problem) {
    map([this, problem](const NodeList::Part &cells) {
        for (auto cell: cells) {
            switch (adaptation_criterion(problem, cell)) {
                case AdaptationFlag::NONE:
                    reset_adaptation_flag(cell);
                    break;

                case AdaptationFlag::SPLIT:
                    set_need_to_split(cell);
                    break;

                default:
                    break;
            }
        }
    });
}

void Mesh::check_flags() {
#ifdef FULL_CHECK
    map([this](const NodeList::Part &cells) {
        for (auto cell: cells) {
            int c_level = (int) cell->level() + static_cast<int>(cell->adaptation_flag());

            // Новый уровень ячейки в пределах нормы
            if (c_level < 0 || c_level > (int)m_max_level) {
                std::cerr << "   cell level: " << cell->level() << "\n";
                std::cerr << "   cell flag: " << static_cast<int>(cell->adaptation_flag()) << "\n";
                throw runtime_error("Cell level out of range");
            }

            // Новые уровни соседей отличаются менее, чем на два
            for (auto &nei: cell->neighbors()) {
                if (m_fake->contain(nei)) {
                    continue;
                }
                int n_level = (int) nei->level() + static_cast<int>(nei->adaptation_flag());

                if (fabs(c_level - n_level) > 1) {
                    NodeList::Ptr bug_cells = NodeList::create(mpi::rank());
                    bug_cells->push_back(cell);
                    bug_cells->push_back(nei);
                    VTKWriter::write(*bug_cells, "check_flags_fail.vtu");
                    std::cerr << "   cell level: " << cell->level() << "\n";
                    std::cerr << "   cell flag: " << static_cast<int>(cell->adaptation_flag()) << "\n";
                    std::cerr << "   neighbor level: " << nei->level() << "\n";
                    std::cerr << "   neighbor flag: " << static_cast<int>(nei->adaptation_flag()) << "\n";
                    throw runtime_error("Adaptation flag imbalance");
                }
            }
        }
    });
#endif
}

void Mesh::set_need_to_coarse(Cell::Ref cell) {
    auto flag = cell->level() > 0 ? AdaptationFlag::COARSE : AdaptationFlag::NONE;
    cell->set_adaptation_flag(flag);
}

void Mesh::reset_adaptation_flag(Cell* cell) {
    if (cell->adaptation_flag() != AdaptationFlag::COARSE)
        return;

    // мы здесь, а значит ячейка собиралась огрубиться,

    // необходимо оповестить её сиблингов, что ничего не выйдет
    cell->set_adaptation_flag(AdaptationFlag::NONE);
    for(auto& sibling: cell->siblings()) {
        reset_adaptation_flag(sibling);
    }

    // а теперь оповестим соседей, чтобы они сильно не отставали
    for(auto& nei: cell->neighbors()) {
        // сосед имеет уровень ниже и хочет огрубиться,
        // во всех остальных случаях баланс 1:2 сохранится
        if (nei->adaptation_flag() == AdaptationFlag::COARSE &&
            nei->level() < cell->level()) {
            reset_adaptation_flag(nei);
        }
    }
}

void Mesh::reset_adaptation_flag(Cell::Ref cell) {
    reset_adaptation_flag(cell.get());
}

void Mesh::set_need_to_split(Cell::Ref cell) {
    set_need_to_split(cell.get());
}

void Mesh::set_need_to_split(Cell* cell) {
    if (cell->adaptation_flag() == AdaptationFlag::SPLIT)
        return;

    // изначально ячейка хотела огрубиться
    if (cell->adaptation_flag() == AdaptationFlag::COARSE) {
        cell->set_adaptation_flag(AdaptationFlag::SPLIT);

        // надо оповестить сиблингов, что ничего не выйдет
        for(auto& sib: cell->siblings()) {
            reset_adaptation_flag(sib);
        }
    }
    cell->set_adaptation_flag(AdaptationFlag::SPLIT);


    // теперь оповестим соседей, чтобы сильно не отставали
    for (auto &nei: cell->neighbors()) {
        if (nei->level() + static_cast<int>(nei->adaptation_flag()) < cell->level()) {
            // сосед после адаптации будет отставать

            if (nei->level() < cell->level()) {
                // сосед ниже уровнем
                set_need_to_split(nei);
            } else {
                // сосед такого же уровня
                reset_adaptation_flag(nei);
            }
        }
    }
}

void Mesh::balance_flags() {
    for (auto &ghost: m_ghosts) {
        for (auto cell: *ghost->cells()) {
            auto recv_flag = cell->adaptation_flag();
            // надо сбросить флаг, иначе автобалансировка не сработает
            set_need_to_coarse(cell);

            switch (recv_flag) {
                case AdaptationFlag::COARSE:
                    set_need_to_coarse(cell);
                    break;

                case AdaptationFlag::NONE:
                    reset_adaptation_flag(cell);
                    break;

                case AdaptationFlag::SPLIT:
                    set_need_to_split(cell);
                    break;
            }
        }
    }
}

void Mesh::apply_flags(Problem* problem, bool debug) {
    vector<vector<Cell::Ptr>> cells_to_refine(m_max_level + 1);
    vector<vector<Cell::Ptr>> cells_to_coarse(m_max_level + 1);

    // В этом блоке резервируем места для ячеек
    // Это не обязательно, но может увеличить производительность
    {
        vector<size_t> refine_count(m_max_level + 1, 0);
        vector<size_t> coarse_count(m_max_level + 1, 0);

        for(auto cell: *m_cells) {
            switch (cell->adaptation_flag()) {
                case AdaptationFlag::NONE:
                    break;
                case AdaptationFlag::SPLIT:
                    ++refine_count[cell->level()];
                    break;
                default:
                    ++coarse_count[cell->level()];
                    break;
            }
        }
        for(auto &ghost: m_ghosts) {
            for(auto cell: *ghost->cells()) {
                switch (cell->adaptation_flag()) {
                    case AdaptationFlag::NONE:
                        break;
                    case AdaptationFlag::SPLIT:
                        ++refine_count[cell->level()];
                        break;
                    default:
                        ++coarse_count[cell->level()];
                        break;
                }
            }
        }
        for(size_t level = 0; level <= m_max_level; ++level) {
            cells_to_refine[level].reserve(refine_count[level]);
            cells_to_coarse[level].reserve(coarse_count[level]);
        }
    }


    for(auto cell: *m_cells) {
        switch (cell->adaptation_flag()) {
            case AdaptationFlag::NONE:
                break;
            case AdaptationFlag::SPLIT:
                cells_to_refine[cell->level()].emplace_back(cell);
                break;
            default:
                cells_to_coarse[cell->level()].emplace_back(cell);
                break;
        }
    }

    for(auto &ghost: m_ghosts) {
        for(auto cell: *ghost->cells()) {
            switch (cell->adaptation_flag()) {
                case AdaptationFlag::NONE:
                    break;
                case AdaptationFlag::SPLIT:
                    cells_to_refine[cell->level()].emplace_back(cell);
                    break;
                default:
                    cells_to_coarse[cell->level()].emplace_back(cell);
                    break;
            }
        }
    }


    // Проходим с нижнего уровня до верхнего и разбиваем ячейки
    for (size_t level = 0; level < m_max_level; ++level) {
        for (auto& cell: cells_to_refine[level]) {
            split_cell(cell, false, problem);
        }
    }
    for (size_t level = m_max_level; level > 0; --level) {
        for (auto &cell: cells_to_coarse[level]) {
            if (cell->parent() && cell->parent()->can_coarse()) {
                coarse_cell(cell->parent(), problem);
            }
        }
    }

    for(auto it = m_fake->begin(); it != m_fake->end(); ) {
        auto cell = *it;
        ++it;
        if (m_fake->contain(cell) && !cell->is_childless()) {
            m_fake->erase(cell);
        }
    }
}

#if DIM2

void Mesh::split_cell(Cell::Ref p_cell, bool ignore_neighbors, Problem* problem) {
#ifdef FULL_CHECK
    // Предполагается, что разбиение ячеек проводится последовательно,
    // начиная с нулевого уровня адаптации и до максимального,
    // поэтому необходимость в следующем коде отпадает

    if (!ignore_neighbors) {
        // Баланс уровней 2:1 между соседями сохраняется на каждом этапе.
        // Поэтому разбиваем сначала соседей меньшего уровня.
        for (auto &face: p_cell->faces()) {
            if (!face->neighbor_exist(p_cell)) {
                continue;
            }

            auto neighbor = face->neighbor(p_cell);
            auto flag = face->flag();

            if (flag == FaceFlag::ORDER && p_cell->level() > neighbor->level()) {
                // Ранее данная строка была. Зачем?
                if (neighbor->is_remote()) {
                    continue;
                }

                std::cerr << "Priority of refinement is not respected\n";
                split_cell(neighbor, false, problem);
            }
        }
    }
#endif

    // Я пробовал делать всё в индексной нотации, но это не упрощает запись,
    // поэтому в дальнейшем для именования вершин и граней используются
    // сокращения l(eft), r(ight), b(ottom), t(op).

    // Создаем вершины
    auto face_center = [](Cell::Ref cell, Side side) -> Vertex::Ptr {
        auto faces = cell->faces(side);
        return faces->is_complex() ? faces->at(0)->vertex(1) : Vertex::create(faces->center());
    };

    // lbv - left-bottom vertex, ccv - center-center vertex
    Vertex::Ptr lbv = p_cell->corner(Side::LEFT, Side::BOTTOM);
    Vertex::Ptr lcv = face_center(p_cell, Side::LEFT);
    Vertex::Ptr ltv = p_cell->corner(Side::LEFT, Side::TOP);

    Vertex::Ptr cbv = face_center(p_cell, Side::BOTTOM);
    Vertex::Ptr ccv = Vertex::create_2d(p_cell->center_2d());
    Vertex::Ptr ctv = face_center(p_cell, Side::TOP);

    Vertex::Ptr rbv = p_cell->corner(Side::RIGHT, Side::BOTTOM);
    Vertex::Ptr rcv = face_center(p_cell, Side::RIGHT);
    Vertex::Ptr rtv = p_cell->corner(Side::RIGHT, Side::TOP);


    // Создаем грани
    // Вертикальные грани
    Face::Ptr lbf = Face::create(lbv, lcv);
    Face::Ptr ltf = Face::create(lcv, ltv);

    Face::Ptr cbf = Face::create(cbv, ccv);
    Face::Ptr ctf = Face::create(ccv, ctv);

    Face::Ptr rbf = Face::create(rbv, rcv);
    Face::Ptr rtf = Face::create(rcv, rtv);

    // Горизонтальные грани
    Face::Ptr brf = Face::create(rbv, cbv);
    Face::Ptr blf = Face::create(cbv, lbv);

    Face::Ptr crf = Face::create(rcv, ccv);
    Face::Ptr clf = Face::create(ccv, lcv);

    Face::Ptr trf = Face::create(rtv, ctv);
    Face::Ptr tlf = Face::create(ctv, ltv);

    // Внешние грани по направлению должны совпадать с родительскими
    if (p_cell->faces(Side::LEFT)->outward_normal().dot(lbf->outward_normal()) < 0.0) {
        lbf->reverse();
        ltf->reverse();
    }
    if (p_cell->faces(Side::BOTTOM)->outward_normal().dot(blf->outward_normal()) < 0.0) {
        blf->reverse();
        brf->reverse();
    }
    if (p_cell->faces(Side::RIGHT)->outward_normal().dot(rbf->outward_normal()) < 0.0) {
        rbf->reverse();
        rtf->reverse();
    }
    if (p_cell->faces(Side::TOP)->outward_normal().dot(tlf->outward_normal()) < 0.0) {
        tlf->reverse();
        trf->reverse();
    }

    // Создадим ячейки на этих гранях
    vector<Cell::Ptr> cells = {
            Cell::create(lbf, blf, cbf, clf),
            Cell::create(cbf, brf, rbf, crf),
            Cell::create(ltf, clf, ctf, tlf),
            Cell::create(ctf, crf, rtf, trf)
    };

    for (size_t k = 0; k < 4; ++k) {
        cells[k]->set_level(p_cell->level() + 1);
        cells[k]->set_parent(p_cell);
        cells[k]->set_owner(p_cell->owner());

        cells[k]->set_z(4 * p_cell->z() + k);
        cells[k]->set_id(z_curve_offset_2d(p_cell->level()) + cells[k]->z());

        if (!problem) {
            cells[k]->copy_data_from(p_cell);
        }
    }


    // Перенос данных
    if (problem) {
        problem->split_data(p_cell, cells);
    }

    // debug
    if (!p_cell->is_childless()) {
        // Данная ситуация может возникать когда пришедшие ячейки уже были у нас, но
        // в какой то момент были переданы другому процессу. В идеале необходимо
        // удалять ячейки и из иерархии тоже (не только из списков) но пока что
        // опробуем более простую стратегию.
        p_cell->remove_children();
    }
    // debug

    // Установим детей
    p_cell->set_children(cells);

    // Свяжем сиблингов
    if (p_cell->is_local()) {
        ctf->set_inward_cell(cells[2]);
        ctf->set_outward_cell(cells[3]);

        cbf->set_inward_cell(cells[0]);
        cbf->set_outward_cell(cells[1]);

        clf->set_inward_cell(cells[0]);
        clf->set_outward_cell(cells[2]);

        crf->set_inward_cell(cells[1]);
        crf->set_outward_cell(cells[3]);
    }

    // Свяжем дочерние ячейки с соседями
    if (!ignore_neighbors) {
        connect_children_with_neighbors(p_cell);
        children_to_lists(p_cell);
    }

    p_cell->set_adaptation_flag(AdaptationFlag::NONE);
}

#else

template <class T>
void rotate(T& t1, T& t2, T& t3, T& t4, int r) {
    int R = (r + 400) % 4;
    if (R < 2) {
        if (R > 0) {
            // R == 1
            T temp = std::move(t4);
            t4 = std::move(t3);
            t3 = std::move(t2);
            t2 = std::move(t1);
            t1 = std::move(temp);
            return;
        }
    }
    else {
        if (R < 3) {
            // R == 2
            T temp = std::move(t1);
            t1 = std::move(t3);
            t3 = std::move(temp);
            temp = std::move(t2);
            t2 = std::move(t4);
            t4 = std::move(temp);
        }
        else {
            // R == 3
            T temp = std::move(t1);
            t1 = std::move(t2);
            t2 = std::move(t3);
            t3 = std::move(t4);
            t4 = std::move(temp);
            return;
        }
    }
}

void Mesh::split_cell(Cell::Ref p_cell, bool ignore_neighbors, Problem* problem) {
    auto get_side = [p_cell](int dir, int c) -> Side {
        switch (dir) {
            case 0:
                return (c < 0 ? Side::LEFT : Side::RIGHT);
            case 1:
                return (c < 0 ? Side::BOTTOM : Side::TOP);
            default:
                return (c < 0 ? Side::BACK : Side::FRONT);
        }
    };

    auto get_faces = [p_cell, get_side](int dir, int c) -> FacesList::Ptr {
        return p_cell->faces(get_side(dir, c));
    };

    // Введем локальную систему координат
    Vector3d center = p_cell->center();
    array<Vector3d, 3> axis;
    for (int dir = 0; dir < 3; ++dir) {
        axis[dir] = get_faces(dir, 1)->center() - get_faces(dir, -1)->center();
        axis[dir].normalize();
    }

    array<array<array<Vertex::Ptr, 3>, 3>, 3> v_grid;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                v_grid[i][j][k] = nullptr;
            }
        }
    }

    array<array<array<Face::Ptr, 2>, 2>, 9> f_grid;
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                f_grid[i][j][k] = nullptr;
            }
        }
    }

    for (int dir: {0, 1, 2}) {
        for (int c: {-1, 1}) {
            auto faces = get_faces(dir, c);

            array<array<Vertex::Ptr, 3>, 3> grid;
            array<array<Face::Ptr, 2>, 2> sub_faces;

            if (faces->is_simple()) {
                grid[0][0] = faces->at(0)->vertex(0);
                grid[2][0] = faces->at(0)->vertex(1);
                grid[2][2] = faces->at(0)->vertex(2);
                grid[0][2] = faces->at(0)->vertex(3);

                grid[0][1] = Vertex::create(0.5 * (grid[0][0]->v() + grid[0][2]->v()));
                grid[1][0] = Vertex::create(0.5 * (grid[0][0]->v() + grid[2][0]->v()));
                grid[1][2] = Vertex::create(0.5 * (grid[0][2]->v() + grid[2][2]->v()));
                grid[2][1] = Vertex::create(0.5 * (grid[2][0]->v() + grid[2][2]->v()));

                grid[1][1] = Vertex::create(faces->center());

                sub_faces[0][0] = Face::create(grid[0][0], grid[1][0], grid[1][1], grid[0][1]);
                sub_faces[1][0] = Face::create(grid[1][0], grid[2][0], grid[2][1], grid[1][1]);
                sub_faces[1][1] = Face::create(grid[1][1], grid[2][1], grid[2][2], grid[1][2]);
                sub_faces[0][1] = Face::create(grid[0][1], grid[1][1], grid[1][2], grid[0][2]);
            } else {
                grid[0][0] = faces->at(0)->vertex(0);
                grid[2][0] = faces->at(1)->vertex(1);
                grid[2][2] = faces->at(2)->vertex(2);
                grid[0][2] = faces->at(3)->vertex(3);

                grid[0][1] = faces->at(0)->vertex(3);
                grid[1][0] = faces->at(0)->vertex(1);
                grid[1][2] = faces->at(2)->vertex(3);
                grid[2][1] = faces->at(1)->vertex(2);

                grid[1][1] = faces->at(0)->vertex(2);

                sub_faces[0][0] = faces->at(0);
                sub_faces[1][0] = faces->at(1);
                sub_faces[1][1] = faces->at(2);
                sub_faces[0][1] = faces->at(3);
            }

            // индекс правильной начальной вершины
            int zero_idx = 0;
            for (int i = 1; i < 4; ++i) {
                // Угловая вершина в локальной системе координат
                Vector3d coord = (faces->is_simple() ? faces->at(0)->vertex(i)->v() :
                                  faces->at((uint) i)->vertex(i)->v()) - center;
                coord = {coord.dot(axis[0]), coord.dot(axis[1]), coord.dot(axis[2])};

                if (coord[(dir + 1) % 3] < 0.0 && coord[(dir + 2) % 3] < 0.0) {
                    zero_idx = i;
                    break;
                }
            }

            // Поворачиваем вершины и грани
            rotate(grid[0][0], grid[2][0], grid[2][2], grid[0][2], -zero_idx);
            rotate(grid[1][0], grid[2][1], grid[1][2], grid[0][1], -zero_idx);
            rotate(sub_faces[0][0], sub_faces[1][0], sub_faces[1][1], sub_faces[0][1], -zero_idx);

            // Реверс вершин и граней, если внешняя нормаль не совпадает с направлением
            if (faces->outward_normal().dot(axis[dir]) < 0.0) {
                std::swap(grid[0][1], grid[1][0]);
                std::swap(grid[2][0], grid[0][2]);
                std::swap(grid[2][1], grid[1][2]);

                std::swap(sub_faces[1][0], sub_faces[0][1]);
            }

            // Заполняем вершины и грани
            for (int ii = 0; ii < 3; ++ii) {
                for (int jj = 0; jj < 3; ++jj) {
                    int i = c + 1;
                    int j = i;
                    int k = i;

                    switch (dir) {
                        case 0:
                            j = ii;
                            k = jj;
                            break;
                        case 1:
                            i = jj;
                            k = ii;
                            break;
                        default:
                            i = jj;
                            j = ii;
                            break;
                    }

                    if (!v_grid[i][j][k]) {
                        v_grid[i][j][k] = grid[ii][jj];
                    }
                }
            }
            for (int ii = 0; ii < 2; ++ii) {
                for (int jj = 0; jj < 2; ++jj) {
                    int idx = 3 * dir + c + 1;
                    int i(ii), j(jj);
                    if (dir == 1) {
                        i = jj;
                        j = ii;
                    }

                    if (!f_grid[idx][i][j]) {
                        f_grid[idx][i][j] = sub_faces[ii][jj];
                    }
                }
            }
        }
    }

    v_grid[1][1][1] = Vertex::create(0.5 * (v_grid[0][1][1]->v() + v_grid[2][1][1]->v()));

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            f_grid[1][i][j] = Face::create(v_grid[1][i][j],
                                           v_grid[1][i][j + 1],
                                           v_grid[1][i + 1][j + 1],
                                           v_grid[1][i + 1][j]);

            f_grid[4][i][j] = Face::create(v_grid[i][1][j],
                                           v_grid[i][1][j + 1],
                                           v_grid[i + 1][1][j + 1],
                                           v_grid[i + 1][1][j]);

            f_grid[7][i][j] = Face::create(v_grid[i][j][1],
                                           v_grid[i + 1][j][1],
                                           v_grid[i + 1][j + 1][1],
                                           v_grid[i][j + 1][1]);
        }
    }

    vector<Cell::Ptr> cells;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                auto cell = Cell::create(
                        f_grid[0 + i][j][k],
                        f_grid[3 + j][i][k],
                        f_grid[1 + i][j][k],
                        f_grid[4 + j][i][k],
                        f_grid[7 + k][i][j],
                        f_grid[6 + k][i][j]
                );
                cells.push_back(cell);
            }
        }
    }

    for (size_t k = 0; k < 8; ++k) {
        cells[k]->set_level(p_cell->level() + 1);
        cells[k]->set_parent(p_cell);
        cells[k]->set_owner(p_cell->owner());

        cells[k]->set_z(8 * p_cell->z() + k);
        cells[k]->set_id(z_curve_offset_3d(p_cell->level()) + cells[k]->z());

        if (!problem) {
            cells[k]->copy_data_from(p_cell);
        }
    }


    // Перенос данных
    if (problem) {
        problem->split_data(p_cell, cells);
    }

    // debug
    if (!p_cell->is_childless()) {
        // Данная ситуация может возникать когда пришедшие ячейки уже были у нас, но
        // в какой то момент были переданы другому процессу. В идеале необходимо
        // удалять ячейки и из иерархии тоже (не только из списков) но пока что
        // опробуем более простую стратегию.
        p_cell->remove_children();
    }
    // debug

    // Установим детей
    p_cell->set_children(cells);

    // Свяжем сиблингов
    if (p_cell->is_local()) {
        for (auto &cell: cells) {
            for (auto &face: cell->faces()) {
                if ((face->center() - cell->center()).dot(face->outward_normal()) > 0.0) {
                    face->set_inward_cell(cell);
                } else {
                    face->set_outward_cell(cell);
                }
            }
        }
    }

    // Свяжем дочерние ячейки с соседями
    if (!ignore_neighbors) {
        connect_children_with_neighbors(p_cell);
        children_to_lists(p_cell);
    }

    p_cell->set_adaptation_flag(AdaptationFlag::NONE);
}

#endif

void Mesh::connect_children_with_neighbors(Cell::Ref cell) {
#if DIM2
    auto side_list = TWO_DIMENSIONAL_SIDES;
#else
    auto side_list = THREE_DIMENSIONAL_SIDES;
#endif

    for(auto side: side_list) {
        auto faces = cell->faces(side);

        // Ячейка на границе
        if (faces->single()->flag() != FaceFlag::ORDER) {
            connect_children_with_fake_cells(cell, side);
        } else {
            // Обычная внутренняя ячейка
            if (faces->is_simple()) {
                connect_children_with_single_cell(cell, side);
            } else {
                connect_children_with_two_cells(cell, side);
            }
        }
    }
}

void Mesh::connect_children_with_fake_cells(Cell::Ref cell, Side side) {
    // Для гост ячеек граница не хранится
    if (cell->is_remote()) {
        return;
    }

    // На границе должна быть единственная ячейка,
    // фейковые ячейки вперед сетки не адаптируются
    auto old_neighbor = cell->faces(side)->single()->neighbor(cell);
    auto children = cell->children_by_side(side);

    static_vector<Face::Ptr, 4> faces;
    static_vector<Cell::Ptr, 4> neighbors;
    for(auto& child: children) {
        faces.push_back(child->faces(side)->single());
        if (child->is_correct_face(side)) {
            faces.back()->set_inward_cell(child);
        }
        else {
            faces.back()->set_outward_cell(child);
        }
        neighbors.push_back(create_neighbor(child, side));
    }

    for(uint i = 0; i < children.size(); ++i) {
        neighbors[i]->set_level(children[i]->level());
        neighbors[i]->set_faces(opposite_side(side), FacesList::create(faces[i]));
        connect_cells(side, children[i], neighbors[i]);

        FaceFlag flag = cell->faces(side)->single()->flag();
        faces[i]->set_flag(flag);
        if (flag == FaceFlag::PERIODIC) {
            neighbors[i]->set_data(old_neighbor->data_holder());
        }
        else if (flag == FaceFlag::ZOE) {
            neighbors[i]->set_data(children[i]->data_holder());
        }
        else {
            neighbors[i]->copy_data_from(old_neighbor);
        }

        m_fake->push_back(neighbors[i]);
    }

    m_fake->erase(old_neighbor);
}

void Mesh::connect_children_with_single_cell(Cell::Ref cell, Side side) {
    auto children = cell->children_by_side(side);
    static_vector<Face::Ptr, 4> faces;
    for (auto &child: children) {
        faces.push_back(child->faces(side)->single());
    }

    auto face = cell->faces(side)->single();
    if (!face->neighbor_exist(cell)) {
        return;
    }

    auto neighbor = face->neighbor(cell);
    if (cell->is_remote() && neighbor->is_remote()) {
        face->set_outward_cell(nullptr);
        face->set_inward_cell(nullptr);
        return;
    }
    Side opp_side = cell->opposite_side(side);
    neighbor->set_faces(opp_side, FacesList::create(faces));

    for (auto &child: children) {
        connect_cells(side, child, neighbor);
    }
}

void Mesh::connect_children_with_two_cells(Cell::Ref cell, Side side) {
    auto children = cell->children_by_side(side);

    for(uint i = 0; i < cell->faces(side)->size(); ++i) {
        auto face = cell->faces(side)->at(i);

        if (!face->neighbor_exist(cell)) {
            continue;
        }
        auto neighbor = face->neighbor(cell);
        if (cell->is_remote() && neighbor->is_remote()) {
            face->set_outward_cell(nullptr);
            face->set_inward_cell(nullptr);
            continue;
        }
        Side opp_side = cell->opposite_side(side);
        auto opp_face = neighbor->faces(opp_side)->single();

        double eps = 1.0e-8 * (opp_face->vertex(1)->v() - opp_face->vertex(0)->v()).squaredNorm();

        for (uint idx = 0; idx < children.size(); ++idx) {
            if ((opp_face->center() - children[idx]->faces(side)->center()).squaredNorm() < eps) {
                neighbor->set_faces(opp_side, children[idx]->faces(side));
                connect_cells(side, children[idx], neighbor);
                break;
            }
        }
    }
}

void Mesh::children_to_lists(Cell::Ref cell) {
    for(auto list: cell->get_lists()) {

        switch (list->role()) {
            case ListRole::MESH:
                // Если родитель находится в основной сетке, то и все дети тоже
                for(auto& child: cell->children()) {
                    list->push_forward(child);
                }
                break;

            case ListRole::BORDER:
                // Если данный ребёнок в граничном слое а сосед удалённый   - добавляем ребёнка
                for(auto& child: cell->children()) {
                    if (list->contain(child)) { continue; }

                    for (auto &neighbor: child->neighbors()) {
                        if (m_fake->contain(neighbor)) { continue; }

                        if (neighbor->has_role(ListRole::GHOST)) {
                            if (neighbor->owner() == list->rank()) {
                                list->push_forward(child);
                            }
                        }
                    }
                }
                break;

            case ListRole::GHOST:
                // Если данный ребёнок в гост слоё а сосед локальный - добавляем ребёнка
                for(auto& child: cell->children()) {
                    if (list->contain(child)) { continue; }

                    for (auto &neighbor: child->neighbors()) {
                        if (m_fake->contain(neighbor)) { continue; }

                        if (neighbor->has_role(ListRole::MESH)) {
                            list->push_forward(child);
                        }
                    }
                }
                break;
        }
    }
    cell->remove_from_all_lists();
}

void Mesh::coarse_cell(Cell::Ref cell, Problem* problem) {
#if DIM2
    auto side_list = TWO_DIMENSIONAL_SIDES;
#else
    auto side_list = THREE_DIMENSIONAL_SIDES;
#endif

    // Следующие две ситуации возникают при огрублении гост слоя. Часто
    // после балансировки у гост слоя сохраняется иерархия, поэтому
    // формально огрубление не может произойти.
    // Проще говоря, есть супер-гост слой ячеек, который недоступен из
    // каких-либо списков, добраться до него можно только по указателям
    // из обычного гост-слоя, так вот для этих супер-гост ячеек также
    // необходимо поддерживать баланс 1:2, иначе может возникнуть баг
    // при огрублении
    for (auto &child: cell->children()) {
        if (!child->is_childless()) {
            child->remove_children();
        }

        for (auto side: side_list) {
            auto faces = child->faces(side);

            if (!faces->is_simple()) {
#if DIM2
                Face::Ptr single_face = Face::create(faces->vertex(0), faces->vertex(1));
#else
                Face::Ptr single_face = Face::create(faces->vertex(0), faces->vertex(1),
                                                     faces->vertex(2), faces->vertex(3));
#endif

                single_face->set_flag(faces->single()->flag());
                if (child == faces->single()->inward_cell()) {
                    single_face->set_inward_cell(child);
                } else {
                    single_face->set_outward_cell(child);
                }
                child->set_faces(side, FacesList::create(single_face));
            }
        }
    }

    // Запишем данные
    size_t cell_data_size = 0;
    if (problem) {
        cell_data_size = problem->cell_data_size();
    } else {
        for (auto &some_child: cell->children()) {
            size_t some_cell_data_size = some_child->data_holder()->size();
            if (some_child != nullptr &&
                some_cell_data_size > 0) {
                cell_data_size = some_cell_data_size;
                break;
            }
        }
    }
    cell->data_holder()->resize(cell_data_size);
    if (!cell->is_remote()) {
        problem->coarse_data(cell, cell->children());
    }

    // Создадим стороны и связи
    for (auto side: side_list) {
        auto children = cell->children_by_side(side);

        // Если бы у ячейки было более одного соседа, то мы бы сюда не попали,
        // так что берём последнюю грань (она одна)
        static_vector<Face::Ptr, 4> faces;
        for (auto &child: children) {
            faces.push_back(child->faces(side)->single());
        }

        // В зависимости от стороны обновим ссылку на грани
        for (auto &face: faces) {
            if (cell->is_correct_face(face)) {
                face->set_inward_cell(cell);
            } else {
                face->set_outward_cell(cell);
            }
        }

        // Соседи
        static_vector<Cell::Ptr, 4> neis;    // Все соседи
        Cell::Ptr neighbor = nullptr;        // Существует ненулевой сосед
        bool null_nei = false;               // Существует нулевой сосед

        for (auto &face: faces) {
            if (face->neighbor_exist(cell)) {
                neis.push_back(face->neighbor(cell));
                neighbor = face->neighbor(cell);
            } else {
                neis.push_back(nullptr);
                null_nei = true;
            }
        }

        // Соседи полностью отсутствуют, такое возможно для ячеек
        // на далекой границе гост слоя
        if (!neighbor) {
            continue;
        }

        // соседи есть не у всех граней, скорее всего гост слой
        // с неполноценной ячейкой без внешних связей,
        // поэтому добавим недостающие связи
        // Похоже на костыль, но у Жоры тут всё на костылях, не?
        if (null_nei) {
            if (neighbor->level() == cell->level()) {
                if (cell->is_correct_face(faces[0])) {
                    for (size_t i = 0; i < faces.size(); ++i) {
                        if (!neis[i]) {
                            faces[i]->set_outward_cell(neighbor);
                        }
                    }
                } else {
                    for (size_t i = 0; i < faces.size(); ++i) {
                        if (!neis[i]) {
                            faces[i]->set_inward_cell(neighbor);
                        }
                    }
                }
            } else {
                continue;
            }
        }

        // Проверим, является ли сосед справа один для обоих детей
        if (neis[0] == neis[1]) {
            // Сосед является одним для обоих детей

            // Функция create гарантирует правильное расположение подграней
            auto complex_face = FacesList::create(faces);

#if DIM2
            Face::Ptr new_face = Face::create(complex_face->at(0)->vertex(0), complex_face->at(1)->vertex(1));
#else
            Face::Ptr new_face = Face::create(complex_face->at(0)->vertex(0), complex_face->at(1)->vertex(1),
                                              complex_face->at(2)->vertex(2), complex_face->at(3)->vertex(3));
#endif

            if (faces[0]->inward_cell() == cell) {
                new_face->set_inward_cell(cell);
                new_face->set_outward_cell(neighbor);
            } else {
                new_face->set_outward_cell(cell);
                new_face->set_inward_cell(neighbor);
            }

            auto new_faces_list = FacesList::create(new_face);

            // Установим новую грань нам и соседу
            cell->set_faces(side, new_faces_list);
            for (auto opp_side: side_list) {
                auto flist = neighbor->faces(opp_side);
                for (size_t i = 0; i < flist->size(); ++i) {
                    if (flist->at(i)->neighbor_exist(neighbor) &&
                        flist->at(i)->neighbor(neighbor) == cell) {
                        neighbor->set_faces(opp_side, new_faces_list);
                    }
                }
            }
        } else {
            // Сосед не является одним для обоих детей
            // Проверим, являются ли они фейковыми
            auto flag = faces[0]->flag();

            if (flag != FaceFlag::ORDER) {
                // Предыдущая грань не удалялась - оставляем все как есть

                Face::Ptr single_face = cell->faces(side)->single();
                if (cell->faces(side)->is_complex()) {
#if DIM2
                    single_face = Face::create(cell->faces(side)->vertex(0),
                                               cell->faces(side)->vertex(1));
#else
                    single_face = Face::create(cell->faces(side)->vertex(0),
                                               cell->faces(side)->vertex(1),
                                               cell->faces(side)->vertex(2),
                                               cell->faces(side)->vertex(3));
#endif
                }

                if (cell->is_correct_face(single_face)) {
                    single_face->set_inward_cell(cell);
                }
                else {
                    single_face->set_outward_cell(cell);
                }

                Cell::Ptr old_neighbor;

                if (single_face->neighbor_exist(cell)) {
                    old_neighbor = single_face->neighbor(cell);
                } else {
                    old_neighbor = create_neighbor(cell, side);
                    connect_cells(side, cell, old_neighbor);

                    single_face->set_flag(flag);

                    m_fake->push_back(old_neighbor);
                }

                if (flag == FaceFlag::WALL) {
                    static_vector<double *, 4> nei_buf;
                    for (auto &nei: neis) {
                        nei_buf.push_back(nei->data_holder()->buffer());
                    }

                    double *nei_buf_old = old_neighbor->data_holder()->buffer();
#if DIM2
                    for (uint i = 0; i < old_neighbor->data_holder()->size(); ++i) {
                        nei_buf_old[i] = 0.5 * (nei_buf[0][i] + nei_buf[1][i]);
                    }
#else
                    for (uint i = 0; i < old_neighbor->data_holder()->size(); ++i) {
                        nei_buf_old[i] = 0.25 * (nei_buf[0][i] + nei_buf[1][i] +
                                                 nei_buf[2][i] + nei_buf[3][i]);
                    }
#endif
                }

                for (auto &nei: neis) {
                    auto list = nei->get_list();
                    list->erase(nei);
                }
            } else {
                // С этой стороны обычные ячейки - добавим себе ссылки на них
                auto new_faces = FacesList::create(faces);
                cell->set_faces(side, new_faces);
            }
        }
    }

    // Удалим детей
    for (auto &child: cell->children()) {
        for (auto &list: child->get_lists()) {
            list->push_forward(cell);
        }
        child->remove_from_all_lists();
    }
    cell->remove_children();

    cell->set_adaptation_flag(AdaptationFlag::NONE);
}

size_t Mesh::level_by_radius(double r) {
    static double r_base = -777.0;
    if (r_base < 0.0) {
        switch (m_geometry->type()) {
            case GeometryType::SECTOR:
                r_base = dynamic_cast<const SectorGenerator *>(m_geometry.get())->r_in();
                break;

            case GeometryType::REAL_SECTOR:
                r_base = 2 * dynamic_cast<const RealSectorGenerator *>(m_geometry.get())->r_min();
                break;

            default:
                runtime_error("Cascade adaptation for unknown geometry");
        }
    }

    bool simple_cascade = false;
    size_t level = 0;

    if (r < r_base) {
        return 0;
    }

    if (simple_cascade) {
        // Необходимый уровень разрешения level удовлетворяет условиям
        // 2^(level - 1) < r / r_in < 2^level)
        level = static_cast<size_t>(floor(log2(r / r_base))); // >= 0
    }
    else {
        double r_max = 50.0;
        double n = 2;
        double L0 = m_max_level / pow(log2(r_max / r_base), n);
        level = static_cast<size_t>(floor(L0 * pow(log2(r / r_base), n)));
    }
    return min(level, m_max_level);
}