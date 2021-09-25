//
// Created by 159-mrv on 3/28/18.
//

#include <control/configuration.h>
#include <core/vertex/vertex.h>
#include <core/cell/side.h>
#include <core/cell/cell.h>
#include <core/face/faces_list.h>
#include <core/face/face.h>
#include <core/generator/structured_decomposition.h>
#include <utils/math.h>
#include <core/generator/sector_generator.h>
#include <utils/geom.h>


SectorGenerator::SectorGenerator(const Configuration &config)
        : StructuredGenerator(GeometryType::SECTOR) {

    m_global_r1 = config("mesh", "geometry", "R1").to_double();
    m_global_r2 = config("mesh", "geometry", "R2").to_double();

    m_alpha = config("mesh", "geometry", "Alpha").to_double();

    // В секторе количество ячеек по Y определяется алгоритмически
    m_global_nx = config("mesh", "cells_per_x").to_uint();
    m_global_ny = static_cast<size_t>(math::int_log(m_global_r2 / m_global_r1, 1 + m_alpha / m_global_nx));
    m_global_nz = 1;

    m_type = GeometryType::SECTOR;
}

vector<vector<vector<Vertex::Ptr>>> SectorGenerator::create_vertices(Decomposition *decomp) const {
    double d_alpha = m_alpha / m_global_nx;

    vector<vector<vector<Vertex::Ptr>>> vertices = vector<vector<vector<Vertex::Ptr>>>(m_global_nx + 1);
    for (size_t x = 0; x <= m_global_nx; ++x) {
        vertices[x] = vector<vector<Vertex::Ptr>>(m_global_ny + 1);
        double phi = m_alpha * (double(x) / m_global_nx - 0.5);

        for (size_t y = 0; y <= m_global_ny; ++y) {
            double r;

            if (y == m_global_ny) {
                r = m_global_r2;
            } else if (y == 0) {
                r = m_global_r1;
            } else {
                r = m_global_r1 * pow(1 + d_alpha, y);
            }

            vertices[x][y] = vector<Vertex::Ptr>(1);
            vertices[x][y][0] = Vertex::create(r * sin(phi), r * cos(phi));
        }
    }

    return vertices;
}

Cell::Ptr SectorGenerator::create_base_cell(size_t z) const {
    auto p = u2xy(z);
    size_t x = p[0];
    size_t y = p[1];

    double d_alpha = m_alpha / m_global_nx;

    double phi1 = m_alpha * (double(x    ) / m_global_nx - 0.5);
    double phi2 = m_alpha * (double(x + 1) / m_global_nx - 0.5);

    double r1 = m_global_r1 * pow(1 + d_alpha, y);
    double r2 = m_global_r2;
    if (y < m_global_ny - 1) {
        r2 = m_global_r1 * pow(1 + d_alpha, y + 1);
    }

    auto lt = Vertex::create(r2 * sin(phi1), r2 * cos(phi1));
    auto lb = Vertex::create(r1 * sin(phi1), r1 * cos(phi1));
    auto rt = Vertex::create(r2 * sin(phi2), r2 * cos(phi2));
    auto rb = Vertex::create(r1 * sin(phi2), r1 * cos(phi2));

    auto left_face   = Face::create(lb, lt);
    auto bottom_face = Face::create(rb, lb);
    auto right_face  = Face::create(rb, rt);
    auto top_face    = Face::create(rt, lt);

    auto cell = Cell::create(left_face, bottom_face, right_face, top_face);
    cell->set_z(xyz2u(x, y));
    cell->set_id(cell->z());

    return cell;
}

Face::Ptr SectorGenerator::next_face_by_side(Cell::Ref cell, Side side) const {
    Face::Ptr face_1 = cell->faces(side)->single();
    Vertex::Ptr v1, v2;

    double r1, r2, phi1, phi2;
    if (side == Side::LEFT || side == Side::RIGHT) {
        r1 = face_1->vertex(0)->v().norm();
        r2 = face_1->vertex(1)->v().norm();

        auto cur_c = face_1->center();
        auto opp_c = cell->faces(opposite_side(side))->center();

        phi1 = phi2 = 2 * atan(cur_c[0] / cur_c[1]) - atan(opp_c[0] / opp_c[1]);

    } else {
        phi1 = atan(face_1->vertex(0)->x() / face_1->vertex(0)->y());
        phi2 = atan(face_1->vertex(1)->x() / face_1->vertex(1)->y());

        double cur_r = face_1->vertex(0)->r();
        double opp_r = cell->faces(opposite_side(side))->vertex(0)->r();

        r1 = r2 = 2 * cur_r - opp_r;
    }

    v1 = Vertex::create(r1 * sin(phi1), r1 * cos(phi1));
    v2 = Vertex::create(r2 * sin(phi2), r2 * cos(phi2));

    auto face_2 = Face::create(v1, v2);
    return face_2;
}

void SectorGenerator::print_info(const std::string& tab) const {
    std::cout << tab << "Inner radius: " << std::hash<double>()(m_global_r1) << "\n";
    std::cout << tab << "Outer radius: " << std::hash<double>()(m_global_r2) << "\n";
    std::cout << tab << "Angle:        " << std::hash<double>()(m_alpha) << "\n";
}

double SectorGenerator::r_in() const {
    return m_global_r1;
}

double SectorGenerator::r_out() const {
    return m_global_r2;
}

double SectorGenerator::angle() const {
    return m_alpha;
}