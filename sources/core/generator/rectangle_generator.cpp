//
// Created by 159-mrv on 8/16/18.
//

#include <control/configuration.h>
#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/side.h>
#include <core/cell/cell.h>
#include <core/generator/structured_decomposition.h>
#include <core/face/face.h>
#include <core/generator/rectangle_generator.h>

RectangleGenerator::RectangleGenerator(const Configuration &config)
        : StructuredGenerator(GeometryType::RECTANGLE) {

    m_x_min = config("mesh", "geometry", "x_min").to_double();
    m_y_min = config("mesh", "geometry", "y_min").to_double();

    m_x_max = config("mesh", "geometry", "x_max").to_double();
    m_y_max = config("mesh", "geometry", "y_max").to_double();

    m_x_len = m_x_max - m_x_min;
    m_y_len = m_y_max - m_y_min;

    // В секторе количество ячеек по Y определяется алгоритмически
    m_global_nx = config("mesh", "cells_per_x").to_uint();
    if (config.exist("mesh", "cells_per_y")) {
        m_global_ny = config("mesh", "cells_per_y").to_uint();
    }
    else {
        m_global_ny = static_cast<size_t>(round(m_y_len / (m_x_len / m_global_nx)));
    }
    m_global_nz = 1;

    m_type = GeometryType::RECTANGLE;
}

vector<vector<vector<Vertex::Ptr>>> RectangleGenerator::create_vertices(Decomposition *decomp) const {
    double x1 = m_x_min;
    double y1 = m_y_min;

    double dx = m_x_len / m_global_nx;
    double dy = m_y_len / m_global_ny;

    vector<vector<vector<Vertex::Ptr>>> vertices(m_global_nx + 1);
    for(size_t i = 0; i <= m_global_nx; ++i) {
        vertices[i].resize(m_global_ny + 1);
        for(size_t j = 0; j <= m_global_ny; ++j) {
            vertices[i][j] = vector<Vertex::Ptr>(1);
            vertices[i][j][0] = Vertex::create(x1 + i * dx, y1 + j * dy);
        }
    }

    for(size_t i = 0; i <= m_global_nx; ++i) {
        for(size_t j = 0; j <= m_global_ny; ++j) {
            double r   = vertices[i][j][0]->r();
            double phi = vertices[i][j][0]->phi();
            phi += ROTATION_ANGLE;

            Vector2d v = { r * cos(phi), r * sin(phi) };
            vertices[i][j][0]->move_2d(v);
        }
    }

    return vertices;
}

Cell::Ptr RectangleGenerator::create_base_cell(size_t z) const {
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
    auto bottom_face = Face::create(rb, lb);
    auto right_face  = Face::create(rb, rt);
    auto top_face    = Face::create(rt, lt);

    auto cell = Cell::create(left_face, bottom_face, right_face, top_face);
    cell->set_z(xyz2u(x, y));
    cell->set_id(cell->z());

    return cell;
}

Face::Ptr RectangleGenerator::next_face_by_side(Cell::Ref cell, Side side) const {
    if (cell->faces(side)->is_complex()) {
        throw runtime_error("RectangleGenerator::next_face_by_side(...) error: Complex face");
    }

    auto face = cell->faces(side)->single();

    double h = (face->center() - cell->faces(opposite_side(side))->center()).norm();
    Vector3d n = - h * face->normal(cell);

    Face::Ptr new_face = Face::create(
            Vertex::create(face->vertex(0)->v() + n),
            Vertex::create(face->vertex(1)->v() + n)
    );

    return new_face;
}

void RectangleGenerator::print_info(const std::string& tab) const {
    std::cout << tab << "Global sizes: " << m_global_nx << " x " << m_global_ny << "\n";
    std::cout << tab << "Left-Bottom vertex: (" << x_min() << ", " << y_min() << ")\n";
    std::cout << tab << "Right-Top   vertex: (" << x_max() << ", " << y_max() << ")\n";
    std::cout << tab << "Sizes: " << x_len() << " x " << y_len() << "\n";
}

double RectangleGenerator::x_min() const {
    return m_x_min;
}

double RectangleGenerator::y_min() const {
    return m_y_min;
}

double RectangleGenerator::x_max() const {
    return m_x_max;
}

double RectangleGenerator::y_max() const {
    return m_y_max;
}

double RectangleGenerator::x_len() const {
    return m_x_len;
}

double RectangleGenerator::y_len() const {
    return m_y_len;
}