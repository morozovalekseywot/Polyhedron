//
// Created by 159-mrv on 3/28/18.
//

#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/vertex/vertex.h>
#include <core/cell/side.h>
#include <core/cell/cell.h>
#include <core/generator/structured_decomposition.h>
#include <core/generator/decomposition.h>
#include <utils/math.h>
#include <core/face/face.h>
#include <core/generator/real_sector_generator.h>
#include <utils/geom.h>

RealSectorGenerator::RealSectorGenerator(const Configuration &config)
        : GeometryGenerator(GeometryType::REAL_SECTOR) {

    m_type = GeometryType::REAL_SECTOR;

    m_r1 = config("mesh", "geometry", "R1").to_double();
    m_r2 = config("mesh", "geometry", "R2").to_double();
    m_alpha = config("mesh", "geometry", "Alpha").to_double();

    m_Nx = config("mesh", "cells_per_x").to_uint();
    m_nc = m_alpha < m_limiter ? 1 : 2;
    m_nx = m_Nx / (2 * m_nc);
    if (m_nx == 0) {
        ++m_nx;
    }
    if (2 * m_nc * m_nx != m_Nx) {
        m_Nx = 2 * m_nc * m_nx;
        if (mpi::is_master()) {
            std::cerr << "  Warning: Parameter Nx should be a multiple of the " + to_string(2 * m_nc) + "\n";
            std::cerr << "           Used value Nx = " + to_string(m_Nx) + ".\n";
        }
    }
    m_Ny = static_cast<size_t>(math::int_log(m_r2 / m_r1, 1 + m_alpha / m_Nx));
    m_size = m_nc * m_nx * (m_nx + 2 * m_Ny);
    
    
    create_vertices();
}

void RealSectorGenerator::create_vertices() {
    // Вершины центральной части
    m_center_vertices = vector<vector<Vertex::Ptr>>(m_nc * m_nx + 1);
    if (m_nc < 2) {
        // Вершины базового квадрата
        Vector2d v1 = { 0.0, 0.0 };
        Vector2d v2 = geom::to_cartesian(m_r1, M_PI_2 - 0.5 * m_alpha);
        Vector2d v3 = geom::to_cartesian(m_r1, M_PI_2);
        Vector2d v4 = geom::to_cartesian(m_r1, M_PI_2 + 0.5 * m_alpha);

        // Вектора граней
        Vector2d v0 = v4;
        Vector2d a = v1 - v4;
        Vector2d b = v3 - v4;
        Vector2d c = v2 - v3;
        Vector2d d = v2 - v1;

        Vector2d skew = 0.5 * ((c + d) - (a + b));

        for(size_t i = 0; i <= m_nx; ++i) {
            m_center_vertices[i] = vector<Vertex::Ptr>(m_nx + 1);
            double x = double(i) / m_nx;
            for(size_t j = 0; j <= m_nx; ++j) {
                double y = double(j) / m_nx;
                m_center_vertices[i][j] = Vertex::create_2d(v0 + x * a + y * b + x * y * skew);
            }
        }
    }
    else {
        // Вершины базовых квадратов
        Vector2d v1 = { 0.0, 0.0 };
        Vector2d v2 = geom::to_cartesian(m_r1, M_PI_2);
        Vector2d v3 = geom::to_cartesian(m_r1, M_PI_2 + 0.25 * m_alpha);
        Vector2d v4 = geom::to_cartesian(m_r1, M_PI_2 + 0.50 * m_alpha);
        Vector2d v5 = geom::to_cartesian(m_r1, M_PI_2 - 0.25 * m_alpha);
        Vector2d v6 = geom::to_cartesian(m_r1, M_PI_2 - 0.50 * m_alpha);

        // Вектора первого квадрата
        Vector2d v0 = v4;
        Vector2d a = v1 - v4;
        Vector2d b = v3 - v4;
        Vector2d c = v2 - v3;
        Vector2d d = v2 - v1;

        Vector2d skew = 0.5 * ((c + d) - (a + b));

        for(size_t i = 0; i <= m_nx; ++i) {
            m_center_vertices[i] = vector<Vertex::Ptr>(m_nx + 1);
            double x = double(i) / m_nx;
            for(size_t j = 0; j <= m_nx; ++j) {
                double y = double(j) / m_nx;
                m_center_vertices[i][j] = Vertex::create_2d(v0 + x * a + y * b + x * y * skew);
            }
        }

        // Вектора граней второго квадрата
        v0 = v1;
        a = v6 - v1;
        b = v2 - v1;
        c = v5 - v2;
        d = v5 - v6;

        skew = 0.5 * ((c + d) - (a + b));

        for(size_t i = 0; i <= m_nx; ++i) {
            m_center_vertices[m_nx + i] = vector<Vertex::Ptr>(m_nx + 1);
            double x = double(i) / m_nx;
            for(size_t j = 0; j <= m_nx; ++j) {
                double y = double(j) / m_nx;
                m_center_vertices[m_nx + i][j] = Vertex::create_2d(v0 + x * a + y * b + x * y * skew);
            }
        }
    }

    // Создаем вершины секторальной части
    double d_alpha = m_alpha / m_Nx;
    m_sector_vertices = vector<vector<Vertex::Ptr>>(m_Nx + 1);
    for(size_t i = 0; i <= m_Nx; ++i) {
        m_sector_vertices[i] = vector<Vertex::Ptr>(m_Ny + 1);

        if (i <= m_nx) {
            m_sector_vertices[i][0] = m_center_vertices[0][i];
        }
        else if (i <= 3 * m_nx) {
            m_sector_vertices[i][0] = m_center_vertices[i - m_nx][m_nx];
        }
        else {
            m_sector_vertices[i][0] = m_center_vertices[2 * m_nx][4 * m_nx - i];
        }

        double phi = M_PI_2 - m_alpha * (double(i) / m_Nx - 0.5);
        for(size_t j = 1; j < m_Ny; ++j) {
            double r = m_r1 * pow(1 + d_alpha, j);
            m_sector_vertices[i][j] = Vertex::create_2d(geom::to_cartesian(r, phi));
        }
        m_sector_vertices[i][m_Ny] = Vertex::create_2d(geom::to_cartesian(m_r2, phi));
    }


    // ========================================================================
    //     Сглаживание сетки
    // ========================================================================

    int smooth_iter = 2000;                     // число итераций сглаживания
    double sigma = 0.5;                        // сила сглаживания
    double r_max = min(1.5 * m_r1, 0.5 * m_r2);  // радиус сглаживания секторальной части

    vector<Vertex::Ptr> left = vector<Vertex::Ptr>(m_nx + 1);
    vector<Vertex::Ptr> top = vector<Vertex::Ptr>(m_nc * m_nx + 1);
    vector<Vertex::Ptr> right = vector<Vertex::Ptr>(m_nx + 1);

    for(size_t i = 0; i <= m_Nx; ++i) {
        if (i <= m_nx) {
            left[i] = m_sector_vertices[i][1];
        }
        if (m_nx <= i && i <= 3 * m_nx) {
            top[i - m_nx] = m_sector_vertices[i][1];
        }
        if (3 * m_nx <= i) {
            right[4 * m_nx - i] = m_sector_vertices[i][1];
        }
    }

    size_t j_max = 0;
    while(m_sector_vertices[0][j_max]->r() < r_max) {
        ++j_max;
    }

    for (int k = 0; k < smooth_iter; ++k) {
        Vector2d avg = {0.0, 0.0};

        // Сглаживание секторальной части
        for (size_t j = 1; j < j_max; ++j) {
            double xi = sigma * (1.0 - m_sector_vertices[0][j]->r() / r_max);
            for (size_t i = 1; i < m_Nx; ++i) {
                if (0 < i && i < m_Nx) {
                    avg = (m_sector_vertices[i + 1][j]->v_2d() + m_sector_vertices[i - 1][j]->v_2d() +
                            m_sector_vertices[i][j - 1]->v_2d() + m_sector_vertices[i][j + 1]->v_2d()) / 4.0;
                } else {
                    avg = (m_sector_vertices[i][j - 1]->v_2d() + m_sector_vertices[i][j + 1]->v_2d()) / 2.0;
                }
                Vector2d hmm = (1 - xi) * m_sector_vertices[i][j]->v_2d() + xi * avg;
                m_sector_vertices[i][j]->move_2d(hmm);
            }
        }


        // Сглаживание центральной части
        for (size_t i = 0; i <= m_nc * m_nx + m_nc - 2; ++i) {
            for (size_t j = 1; j <= m_nx; ++j) {
                auto lv = i > 0 ? m_center_vertices[i - 1][j] : left[j];
                auto rv = i < m_nc * m_nx ? m_center_vertices[i + 1][j] : right[j];
                auto bv = m_center_vertices[i][j - 1];
                auto tv = j < m_nx ? m_center_vertices[i][j + 1] : top[i];

                avg = (lv->v_2d() + rv->v_2d() + bv->v_2d() + tv->v_2d()) / 4.0;
                if (tv == lv || tv == rv) {
                    avg = (lv->v_2d() + rv->v_2d() + bv->v_2d()) / 3.0;
                }
                Vector2d hmm = (1 - sigma) * m_center_vertices[i][j]->v_2d() + sigma * avg;
                m_center_vertices[i][j]->move_2d(hmm);
            }
        }
    }
}

Face::Ptr RealSectorGenerator::get_center_hface(size_t x, size_t y) const {
    assert(0 <= x && x < m_nc * m_nx);
    assert(0 <= y && y <= m_nx);

    return Face::create(m_center_vertices[x][y],
                        m_center_vertices[x + 1][y]);
}

Face::Ptr RealSectorGenerator::get_center_vface(size_t x, size_t y) const {
    assert(0 <= x && x <= m_nc * m_nx);
    assert(0 <= y && y <= m_nx);

    if (x <= m_nx) {
        // Для левого основного куба вертикальные грани направлены вверх
        return Face::create(m_center_vertices[x][y],
                            m_center_vertices[x][y + 1]);
    } else {
        // Для правого основного куба -- вниз
        return Face::create(m_center_vertices[x][y + 1],
                            m_center_vertices[x][y]);
    }
}

Face::Ptr RealSectorGenerator::get_sector_hface(size_t x, size_t y) const {
    if (y > 0) {
        return Face::create(m_sector_vertices[x][y], m_sector_vertices[x + 1][y]);
    }
    else {
        if (x < m_nx) {
            return get_center_vface(0, x);
        }
        else if (x < 3 * m_nx) {
            return get_center_hface(x - m_nx, m_nx);
        }
        else {
            return get_center_vface(2 * m_nx, 4 * m_nx - x - 1);
        }
    }
}

Face::Ptr RealSectorGenerator::get_sector_vface(size_t x, size_t y) const {
    return Face::create(m_sector_vertices[x][y], m_sector_vertices[x][y + 1]);
}

vector<Cell::Ptr> RealSectorGenerator::create_cells(Decomposition *decomp) const {
    vector<Cell::Ptr> cells;

    // ========================================================================
    //     СОЗДАЕМ ГРАНИ
    // ========================================================================

    // Создаем горизонтальные грани центральной части
    // Все горизонтальные грани направлены вправо
    vector<vector<Face::Ptr>> hfaces1 =
            vector<vector<Face::Ptr>>(m_nc * m_nx, vector<Face::Ptr>(m_nx + 1));
    for (size_t x = 0; x < m_nc * m_nx; ++x) {
        for (size_t y = 0; y <= m_nx; ++y) {
            hfaces1[x][y] = get_center_hface(x, y);
        }
    }

    // Создаем вертикальные грани центральной части
    vector<vector<Face::Ptr>> vfaces1 =
            vector<vector<Face::Ptr>>(m_nc * m_nx + 1, vector<Face::Ptr>(m_nx));
    // Для правого основного куба вертикальные грани направлены вниз
    for (size_t x = 0; x <= m_nc * m_nx; ++x) {
        for (size_t y = 0; y < m_nx; ++y) {
            vfaces1[x][y] = get_center_vface(x, y);
        }
    }

    // Создаем горизонтальные секторальной части
    // Нижний слой уже создан для центральной части
    vector<vector<Face::Ptr>> hfaces2 =
            vector<vector<Face::Ptr>>(m_Nx, vector<Face::Ptr>(m_Ny + 1));
    for (size_t x = 0; x < m_Nx; ++x) {
        if (x < m_nx) {
            hfaces2[x][0] = vfaces1[0][x];
        }
        else if (x < 3 * m_nx) {
            hfaces2[x][0] = hfaces1[x - m_nx][m_nx];
        }
        else {
            hfaces2[x][0] = vfaces1[2 * m_nx][4 * m_nx - x - 1];
        }

        for (size_t y = 1; y <= m_Ny; ++y) {
            hfaces2[x][y] = get_sector_hface(x, y);
        }
    }

    // Создаем вертикальные грани секторальной части
    vector<vector<Face::Ptr>> vfaces2 =
            vector<vector<Face::Ptr>>(m_Nx + 1, vector<Face::Ptr>(m_Ny));
    for (size_t x = 0; x <= m_Nx; ++x) {
        for (size_t y = 0; y < m_Ny; ++y) {
            vfaces2[x][y] = get_sector_vface(x, y);
        }
    }


    // ========================================================================
    //     СОЗДАЕМ ЯЧЕЙКИ
    // ========================================================================

    auto test_neighbor = [&, vfaces1, vfaces2, hfaces1, hfaces2]
            (const Index& idx, Side side) -> Cell::Ptr {

        Index nei_idx = neighbor(idx, side);

        if (nei_idx.part < 0) {
            return nullptr;
        }

        if (nei_idx.part < 1) {
            auto nl_face = vfaces1[nei_idx.x][nei_idx.y];
            auto nb_face = hfaces1[nei_idx.x][nei_idx.y];
            auto nr_face = vfaces1[nei_idx.x + 1][nei_idx.y];
            auto nt_face = hfaces1[nei_idx.x][nei_idx.y + 1];

            return Cell::create(nl_face, nb_face, nr_face, nt_face);
        } else {
            auto nl_face = vfaces2[nei_idx.x][nei_idx.y];
            auto nb_face = hfaces2[nei_idx.x][nei_idx.y];
            auto nr_face = vfaces2[nei_idx.x + 1][nei_idx.y];
            auto nt_face = hfaces2[nei_idx.x][nei_idx.y + 1];

            return Cell::create(nl_face, nb_face, nr_face, nt_face);
        }
    };

    // Создаем ячейки центральной части
    for (size_t x = 0; x < m_nc * m_nx; ++x) {
        for (size_t y = 0; y < m_nx; ++y) {
            Index idx = {0, x, y};

            auto l_face = vfaces1[x][y];
            auto b_face = hfaces1[x][y];
            auto r_face = vfaces1[x + 1][y];
            auto t_face = hfaces1[x][y + 1];

            auto cell = Cell::create(l_face, b_face, r_face, t_face);
            int owner = decomp->rank(cell->center());

            cell->set_owner(owner);

            bool useless = cell->is_remote();
            if (useless) {
                for(auto side: TWO_DIMENSIONAL_SIDES) {
                    auto test_nei = test_neighbor(idx, side);
                    if (!test_nei) {
                        continue;
                    }

                    if (decomp->rank(test_nei->center()) == mpi::rank()) {
                        useless = false;
                        break;
                    }
                }
            }


            if (useless) {
                continue;
            }

            size_t z = get_z(idx);
            cell->set_z(z);
            cell->set_id(z);

            if (x <= m_nx) {
                l_face->set_outward_cell(cell);
            }
            else {
                l_face->set_inward_cell(cell);
            }
            t_face->set_outward_cell(cell);

            b_face->set_inward_cell(cell);
            if (x < m_nx) {
                r_face->set_inward_cell(cell);
            }
            else {
                r_face->set_outward_cell(cell);
            }

            cells.push_back(cell);
        }
    }


    // Создаем ячейки секторальной части
    for (size_t x = 0; x < m_Nx; ++x) {
        for (size_t y = 0; y < m_Ny; ++y) {
            Index idx = {1, x, y};

            auto l_face = vfaces2[x][y];
            auto b_face = hfaces2[x][y];
            auto r_face = vfaces2[x + 1][y];
            auto t_face = hfaces2[x][y + 1];

            auto cell = Cell::create(l_face, b_face, r_face, t_face);
            int owner = decomp->rank(cell->center());

            cell->set_owner(owner);

            bool useless = cell->is_remote();
            if (useless) {
                for(auto side: TWO_DIMENSIONAL_SIDES) {
                    auto test_nei = test_neighbor(idx, side);
                    if (!test_nei) {
                        continue;
                    }

                    if (decomp->rank(test_nei->center()) == mpi::rank()) {
                        useless = false;
                        break;
                    }
                }
            }

            if (useless) {
                continue;
            }

            size_t z = get_z(idx);
            cell->set_z(z);
            cell->set_id(z);

            l_face->set_outward_cell(cell);
            t_face->set_outward_cell(cell);

            b_face->set_inward_cell(cell);
            r_face->set_inward_cell(cell);

            cells.push_back(cell);
        }
    }

    return cells;
}

Cell::Ptr RealSectorGenerator::create_base_cell(size_t z) const {
    if (z < m_nc * m_nx * m_nx) {
        // Центральная часть
        size_t x = z % (m_nc * m_nx);
        size_t y = z / (m_nc * m_nx);

        auto l_face = get_center_vface(x, y);
        auto b_face = get_center_hface(x, y);
        auto r_face = get_center_vface(x + 1, y);
        auto t_face = get_center_hface(x, y + 1);

        auto cell = Cell::create(l_face, b_face, r_face, t_face);
        cell->set_z(z);
        cell->set_id(z);
        return cell;
    }
    else {
        // Секторальная область
        size_t x = (z - m_nc * m_nx * m_nx) % m_Nx;
        size_t y = (z - m_nc * m_nx * m_nx) / m_Nx;

        auto l_face = get_sector_vface(x, y);
        auto b_face = get_sector_hface(x, y);
        auto r_face = get_sector_vface(x + 1, y);
        auto t_face = get_sector_hface(x, y + 1);

        auto cell = Cell::create(l_face, b_face, r_face, t_face);
        cell->set_z(z);
        cell->set_id(z);
        return cell;
    }
}

Face::Ptr RealSectorGenerator::next_face_by_side(Cell::Ref cell, Side side) const {
    size_t z_base = cell->z_base();

    auto lb = cell->left_bottom_vertex()->v_2d();
    auto lt = cell->left_top_vertex()->v_2d();
    auto rb = cell->right_bottom_vertex()->v_2d();
    auto rt = cell->right_top_vertex()->v_2d();

    Vector2d v1, v2;

    // TODO: Сшивка фейковой границы
    if (z_base < m_nc * m_nx * m_nx) {
        if (m_nc < 2) {
            switch (side) {
                case Side::BOTTOM:
                    v1 = geom::symmetric_point(lt, lb, rb);
                    v2 = geom::symmetric_point(rt, lb, rb);
                    break;
                case Side::RIGHT:
                    v1 = geom::symmetric_point(lb, rb, rt);
                    v2 = geom::symmetric_point(lt, rb, rt);
                    break;
                default:
                    throw runtime_error("Error! Trying to create neighbor by invalid side! #1");
            }
        }
        else {
            switch (side) {
                case Side::BOTTOM:
                    v1 = geom::symmetric_point(lt, lb, rb);
                    v2 = geom::symmetric_point(rt, lb, rb);
                    break;
                default:
                    throw runtime_error("Error! Trying to create neighbor by invalid side! #2");
            }
        }
    }
    else {
        // Секторальная часть
        double r, phi_1, phi_2;
        switch (side) {
            case Side::LEFT:
                v1 = geom::symmetric_point(rb, lb, lt);
                v2 = geom::symmetric_point(rt, lb, lt);
                break;

            case Side::RIGHT:
                v1 = geom::symmetric_point(lb, rb, rt);
                v2 = geom::symmetric_point(lt, rb, rt);
                break;

            case Side::TOP:
                phi_1 = atan2(rt[1], rt[0]);
                phi_2 = atan2(lt[1], lt[0]);
                r = 2 * lt.norm() - lb.norm();

                v1 = geom::to_cartesian(r, phi_2);
                v2 = geom::to_cartesian(r, phi_1);
                break;

            default:
                std::cerr << to_string(side) << " FAIL SIDE\n";
                throw runtime_error("Error! Trying to create real_sector neighbor by invalid side! #1");
        }
    }

    return Face::create(Vertex::create_2d(v1), Vertex::create_2d(v2));
}

size_t RealSectorGenerator::n_cells() const noexcept {
    return m_size;
}


RealSectorGenerator::Index RealSectorGenerator::neighbor(const Index& index, Side side) const {
    Index nei;

    if (index.part < 1) {
        // Центральная область
        nei.part = 0;
        nei.x = index.x;
        nei.y = index.y;

        switch (side) {
            case Side::LEFT:
                if (index.x == 0) {
                    // Попадаем в секторальную часть
                    nei.part = 1;
                    nei.x = index.x;
                    nei.y = 0;
                }
                else {
                    --nei.x;
                }
                return nei;

            case Side::RIGHT:
                if (index.x == 2 * m_nx - 1) {
                    // Попадаем в секторальную часть
                    nei.part = 1;
                    nei.x = m_Nx - index.y - 1;
                    nei.y = 0;
                }
                else {
                    ++nei.x;
                }
                return nei;

            case Side::BOTTOM:
                if (index.y == 0) {
                    nei.part = -1;
                }
                --nei.y;
                return nei;

            case Side::TOP:
                if (index.y == m_nx - 1) {
                    // Попадаем в секторальную часть
                    nei.part = 1;
                    nei.x = index.x + m_nx;
                    nei.y = 0;
                }
                else {
                    ++nei.y;
                }
                return nei;

            default:
                throw runtime_error("RealSectorGenerator: Unknown side");
        }
    }
    else {
        // Секторальная область
        nei.part = 1;
        nei.x = index.x;
        nei.y = index.y;

        switch (side) {
            case Side::LEFT:
                if (index.x > 0) {
                    --nei.x;
                }
                else {
                    nei.part = -1;
                }
                return nei;

            case Side::RIGHT:
                if (index.x < m_Nx - 1) {
                    ++nei.x;
                }
                else {
                    nei.part = -1;
                }
                return nei;

            case Side::TOP:
                if (index.y < m_Ny - 1) {
                    ++nei.y;
                }
                else {
                    nei.part = -1;
                }
                return nei;

            case Side::BOTTOM:
                if (index.y > 0) {
                    --nei.y;
                    return nei;
                }

                // попадаем в центральную область
                nei.part = 0;

                if (index.x < m_nx) {
                    nei.x = 0;
                    nei.y = index.x;
                }
                else if (index.x < 3 * m_nx) {
                    nei.x = index.x - m_nx;
                    nei.y = m_nx - 1;
                }
                else {
                    nei.x = m_nc * m_nx - 1;
                    nei.y = m_nx - (index.x - 3 * m_nx) - 1;
                }
                return nei;

            default:
                throw runtime_error("RealSectorGenerator: Unknown side");
        }
    }

}

size_t RealSectorGenerator::neighbor_z(Cell::Ref cell, Side side) const {
    size_t z = cell->z_base();

    if (z < m_nc * m_nx * m_nx) {
        // Центральная область

        size_t x_base = z % (m_nc * m_nx);
        size_t y_base = z / (m_nc * m_nx);

        switch (side) {
            case Side::LEFT:
                if (x_base == 0) {
                    // Попадаем в секторальную часть
                    x_base = y_base;
                    y_base = 0;
                    return sector_z(x_base, y_base);
                }
                else {
                    --x_base;
                    return center_z(x_base, y_base);
                }

            case Side::RIGHT:
                if (x_base == 2 * m_nx - 1) {
                    // Попадаем в секторальную часть
                    x_base = m_Nx - y_base - 1;
                    y_base = 0;
                    return sector_z(x_base, y_base);
                }
                else {
                    ++x_base;
                    return center_z(x_base, y_base);
                }

            case Side::BOTTOM:
                --y_base;
                return center_z(x_base, y_base);

            case Side::TOP:
                if (y_base == m_nx - 1) {
                    // Попадаем в секторальную часть
                    x_base = x_base + m_nx;
                    y_base = 0;
                    return sector_z(x_base, y_base);
                }
                else {
                    ++y_base;
                    return center_z(x_base, y_base);
                }

            default:
                throw runtime_error("RealSectorGenerator: Unknown side");
        }
    }
    else {
        // Секторальная область
        z -= m_nc * m_nx * m_nx;

        size_t x_base = z % m_Nx;
        size_t y_base = z / m_Nx;

        switch (side) {
            case Side::LEFT:
                --x_base;
                break;

            case Side::RIGHT:
                ++x_base;
                break;

            case Side::TOP:
                ++y_base;
                break;

            case Side::BOTTOM:
                // попадаем в центральную область
                if (x_base < m_nx) {
                    x_base = 0;
                    y_base = x_base;
                }
                else if (x_base < 3 * m_nx) {
                    x_base = x_base - m_nx;
                    y_base = m_nx - 1;
                }
                else {
                    x_base = m_nc * m_nx - 1;
                    y_base = m_nx - (x_base - 3 * m_nx) - 1;
                }
                return center_z(x_base, y_base);

            default:
                throw runtime_error("RealSectorGenerator: Unknown side");
        }

        return sector_z(x_base, y_base);
    }
}

void RealSectorGenerator::print_info(const std::string& tab) const {
    std::cout << tab << "Inner radius: " << std::hash<double>()(m_r1) << "\n";
    std::cout << tab << "Outer radius: " << std::hash<double>()(m_r2) << "\n";
    std::cout << tab << "Angle:        " << std::hash<double>()(m_alpha) << "\n";
}

size_t RealSectorGenerator::center_z(size_t x_base, size_t y_base) const {
    return m_nc * m_nx * y_base + x_base;
}

size_t RealSectorGenerator::sector_z(size_t x_base, size_t y_base) const {
    return m_nc * m_nx * m_nx + y_base * m_Nx + x_base;
}

size_t RealSectorGenerator::get_z(const Index& index) const {
    if (index.part == 0) {
        return center_z(index.x, index.y);
    }
    else {
        return sector_z(index.x, index.y);
    }
}
