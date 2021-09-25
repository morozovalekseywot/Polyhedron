//
// Created by 159-mrv on 8/27/19.
//

#include <control/configuration.h>
#include <core/vertex/vertex.h>
#include <core/cell/side.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/cell.h>
#include <core/generator/structured_decomposition.h>
#include <core/face/face.h>
#include <core/generator/cuboid_generator.h>

#if DIM3

CuboidGenerator::CuboidGenerator(const Configuration &config)
        : StructuredGenerator(GeometryType::RECTANGLE) {

    m_x_min = config("mesh", "geometry", "x_min").to_double();
    m_y_min = config("mesh", "geometry", "y_min").to_double();
    m_z_min = config("mesh", "geometry", "z_min").to_double();

    m_x_max = config("mesh", "geometry", "x_max").to_double();
    m_y_max = config("mesh", "geometry", "y_max").to_double();
    m_z_max = config("mesh", "geometry", "z_max").to_double();

    m_x_len = m_x_max - m_x_min;
    m_y_len = m_y_max - m_y_min;
    m_z_len = m_z_max - m_z_min;

    // В секторе количество ячеек по Y определяется алгоритмически
    m_global_nx = config("mesh", "cells_per_x").to_uint();
    if (config.exist("mesh", "cells_per_y")) {
        m_global_ny = config("mesh", "cells_per_y").to_uint();
    }
    else {
        m_global_ny = static_cast<size_t>(round(m_y_len / (m_x_len / m_global_nx)));
    }

    if (config.exist("mesh", "cells_per_z")) {
        m_global_nz = config("mesh", "cells_per_z").to_uint();
    }
    else {
        m_global_nz = static_cast<size_t>(round(m_z_len / (m_x_len / m_global_nx)));
    }

    m_type = GeometryType::CUBOID;
}

vector<vector<vector<Vertex::Ptr>>> CuboidGenerator::create_vertices(Decomposition *decomp) const {
    double x1 = m_x_min;
    double y1 = m_y_min;
    double z1 = m_z_min;

    double dx = m_x_len / m_global_nx;
    double dy = m_y_len / m_global_ny;
    double dz = m_z_len / m_global_nz;

    auto vertices = vector<vector<vector<Vertex::Ptr>>>(m_global_nx + 1);
    for(size_t i = 0; i <= m_global_nx; ++i) {
        vertices[i] = vector<vector<Vertex::Ptr>>(m_global_ny + 1);
        for(size_t j = 0; j <= m_global_ny; ++j) {
            vertices[i][j] = vector<Vertex::Ptr>(m_global_nz + 1);
            for(size_t k = 0; k <= m_global_nz; ++k) {
                vertices[i][j][k] = Vertex::create(x1 + i * dx, y1 + j * dy, z1 + k * dz);
            }
        }
    }

    return vertices;
}

Cell::Ptr CuboidGenerator::create_base_cell(size_t z) const {
    throw runtime_error("Create Base Cell Cuboid Generator");
    auto p = u2xy(z);
    size_t x = p[0];
    size_t y = p[1];

    double dx = m_x_len / m_global_nx;
    double x1 = m_x_min + x * dx;
    double x2 = x1 + dx;

    double dy = m_y_len / m_global_ny;
    double y1 = m_y_min + y * dy;
    double y2 = y1 + dy;

    auto lt = Vertex::create(x1, y2);
    auto lb = Vertex::create(x1, y1);
    auto rt = Vertex::create(x2, y2);
    auto rb = Vertex::create(x2, y1);

    auto left_face   = Face::create(lb, lt);
    auto bottom_face = Face::create(lb, rb);
    auto right_face  = Face::create(rb, rt);
    auto top_face    = Face::create(lt, rt);

    auto cell = Cell::create(left_face, bottom_face, right_face, top_face);
    cell->set_z(xyz2u(x, y));
    cell->set_id(cell->z());

    return cell;
}

Face::Ptr CuboidGenerator::next_face_by_side(Cell::Ref cell, Side side) const {
    if (cell->faces(side)->is_complex()) {
        throw runtime_error("CuboidGenerator::next_face_by_side(...) error: Complex face");
    }

    auto face = cell->faces(side)->single();

    double h = (face->center() - cell->faces(opposite_side(side))->center()).norm();
    Vector3d n = - h * face->normal(cell);

    Face::Ptr new_face = Face::create(
            Vertex::create(face->vertex(0)->v() + n),
            Vertex::create(face->vertex(1)->v() + n),
            Vertex::create(face->vertex(2)->v() + n),
            Vertex::create(face->vertex(3)->v() + n)
    );

    return new_face;
}

void CuboidGenerator::print_info(const std::string& tab) const {
    std::cout << tab << "Global sizes: " << m_global_nx << " x " << m_global_ny << "\n";
    std::cout << tab << "Left-Bottom vertex: (" << x_min() << ", " << y_min() << ")\n";
    std::cout << tab << "Right-Top   vertex: (" << x_max() << ", " << y_max() << ")\n";
    std::cout << tab << "Sizes: " << x_len() << " x " << y_len() << "\n";
}

double CuboidGenerator::x_min() const {
    return m_x_min;
}

double CuboidGenerator::y_min() const {
    return m_y_min;
}

double CuboidGenerator::z_min() const {
    return m_z_min;
}

double CuboidGenerator::x_max() const {
    return m_x_max;
}

double CuboidGenerator::y_max() const {
    return m_y_max;
}

double CuboidGenerator::z_max() const {
    return m_z_max;
}

double CuboidGenerator::x_len() const {
    return m_x_len;
}

double CuboidGenerator::y_len() const {
    return m_y_len;
}

double CuboidGenerator::z_len() const {
    return m_z_len;
}

#endif