//
// Created by 159-mrv on 8/16/18.
//

#include <control/configuration.h>
#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/cell.h>
#include <core/generator/structured_decomposition.h>
#include <core/generator/sector_generator.h>
#include <core/generator/decomposition.h>
#include <core/generator/rectangle_generator.h>
#include <control/mpi_wrapper.h>

StructuredGenerator::StructuredGenerator(GeometryType type)
    : GeometryGenerator(type) {

}

vector<Cell::Ptr> StructuredGenerator::create_cells(Decomposition *decomp) const {
    auto vertices = create_vertices(decomp);

    bool TWO_DIM = vertices[0][0].size() == 1;

    // Грани с внешней нормалью по оси Ox
    vector<vector<vector<Face::Ptr>>> x_faces(m_global_nx + 1);
    for (size_t x = 0; x <= m_global_nx; ++x) {
        x_faces[x].resize(m_global_ny);
        for (size_t y = 0; y < m_global_ny; ++y) {
            x_faces[x][y].resize(m_global_nz);
#if DIM2
            x_faces[x][y][0] = Face::create(
                    vertices[x][y][0],
                    vertices[x][y + 1][0]
            );
#else
            for (size_t z = 0; z < m_global_nz; ++z) {
                x_faces[x][y][z] = Face::create(
                        vertices[x][y][z],
                        vertices[x][y + 1][z],
                        vertices[x][y + 1][z + 1],
                        vertices[x][y][z + 1]
                );
            }
#endif
        }
    }


    // Грани с внешней нормалью по оси Y
    vector<vector<vector<Face::Ptr>>> y_faces(m_global_nx);
    for (size_t x = 0; x < m_global_nx; ++x) {
        y_faces[x].resize(m_global_ny + 1);
        for (size_t y = 0; y <= m_global_ny; ++y) {
            y_faces[x][y].resize(m_global_nz);
#if DIM2
            y_faces[x][y][0] = Face::create(
                    vertices[x + 1][y][0],
                    vertices[x][y][0]
            );
#else
            for(size_t z = 0; z < m_global_nz; ++z) {
                y_faces[x][y][z] = Face::create(
                        vertices[x][y][z],
                        vertices[x][y][z + 1],
                        vertices[x + 1][y][z + 1],
                        vertices[x + 1][y][z]
                );
            }
#endif
        }
    }


    // Грани с внешней нормалью по оси Z
    vector<vector<vector<Face::Ptr>>> z_faces;
#if DIM3
    z_faces.resize(m_global_nx);
    for (size_t x = 0; x < m_global_nx; ++x) {
        z_faces[x].resize(m_global_ny);
        for (size_t y = 0; y < m_global_ny; ++y) {
            z_faces[x][y].resize(m_global_nz + 1);
            for (size_t z = 0; z <= m_global_nz; ++z) {
                z_faces[x][y][z] = Face::create(
                        vertices[x][y][z],
                        vertices[x + 1][y][z],
                        vertices[x + 1][y + 1][z],
                        vertices[x][y + 1][z]
                );
            }
        }
    }

#endif

    vector<Cell::Ptr> cells;
    for(size_t y = 0; y < m_global_ny; ++y) {
        for (size_t x = 0; x < m_global_nx; ++x) {
            for (size_t z = 0; z < m_global_nz; ++z) {
                array<Face::Ptr, 6> faces = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };

                faces[to_int(Side::LEFT)]   = x_faces[x][y][z];
                faces[to_int(Side::RIGHT)]  = x_faces[x + 1][y][z];
                faces[to_int(Side::BOTTOM)] = y_faces[x][y][z];
                faces[to_int(Side::TOP)]    = y_faces[x][y + 1][z];

                Cell::Ptr cell;

                if (TWO_DIM) {
                    cell = Cell::create(
                            faces[to_int(Side::LEFT)],
                            faces[to_int(Side::BOTTOM)],
                            faces[to_int(Side::RIGHT)],
                            faces[to_int(Side::TOP)]
                    );
                }
                else {
                    faces[to_int(Side::BACK)] = z_faces[x][y][z];
                    faces[to_int(Side::FRONT)]  = z_faces[x][y][z + 1];

                    cell = Cell::create(
                            faces[to_int(Side::LEFT)],
                            faces[to_int(Side::BOTTOM)],
                            faces[to_int(Side::RIGHT)],
                            faces[to_int(Side::TOP)],
                            faces[to_int(Side::FRONT)],
                            faces[to_int(Side::BACK)]
                    );
                }

                cell->set_owner(decomp->rank(cell->center()));

                if (cell->is_remote()) {
                    bool useless = true;
                    for(auto& face: faces) {
                        if (!face) { continue; }

                        Vector3d nei_center = cell->center() + 2 * (face->center() - cell->center());
                        if (decomp->rank(nei_center) == mpi::rank()) {
                            useless = false;
                            break;
                        }
                    }

                    if (useless) {
                        continue;
                    }
                }

                size_t u = xyz2u(x, y, z);
                cell->set_z(u);
                cell->set_id(u);

                faces[to_int(Side::LEFT)]->set_outward_cell(cell);
                faces[to_int(Side::RIGHT)]->set_inward_cell(cell);

                faces[to_int(Side::BOTTOM)]->set_outward_cell(cell);
                faces[to_int(Side::TOP)]->set_inward_cell(cell);

                if (!TWO_DIM) {
                    faces[to_int(Side::BACK)]->set_outward_cell(cell);
                    faces[to_int(Side::FRONT)]->set_inward_cell(cell);
                }

                cells.push_back(cell);
            }
        }
    }

    // Отсутствие ячеек считаем ошибкой
    if (cells.empty()) {
        throw runtime_error("ERROR! Creating mesh without cells");
    }

    return cells;
}

size_t StructuredGenerator::global_nx() const noexcept {
    return m_global_nx;
}

size_t StructuredGenerator::global_ny() const noexcept {
    return m_global_ny;
}

size_t StructuredGenerator::global_nz() const noexcept {
    return m_global_nz;
}

size_t StructuredGenerator::n_cells() const noexcept {
    return m_global_nx * m_global_ny * m_global_nz;
}

size_t StructuredGenerator::neighbor_z(Cell::Ref cell, Side side) const {
    auto p = u2xyz(cell->z_base());

    size_t x_nei = p[0];
    size_t y_nei = p[1];
    size_t z_nei = p[2];
    switch (side) {
        case Side::LEFT:
            --x_nei;
            break;

        case Side::RIGHT:
            ++x_nei;
            break;

        case Side::BOTTOM:
            --y_nei;
            break;

        case Side::TOP:
            ++y_nei;
            break;

        case Side::FRONT:
            --z_nei;
            break;

        case Side::BACK:
            ++z_nei;
            break;

        default:
            throw runtime_error("Unknown side");
    }
    return xyz2u(x_nei, y_nei, z_nei);
}

size_t StructuredGenerator::xyz2u(size_t x, size_t y, size_t z) const {
    return (z * m_global_ny + y) * m_global_nx + x;
}

array<size_t, 2> StructuredGenerator::u2xy(size_t u) const {
    size_t x = u % m_global_nx;
    size_t y = u / m_global_nx;
    return {x, y};
}

array<size_t, 3> StructuredGenerator::u2xyz(size_t u) const {
    size_t x = u % m_global_nx;
    size_t zy = u / m_global_nx;
    size_t y = zy % m_global_ny;
    size_t z = zy / m_global_ny;

    return {x, y, z};
}