#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/vertex/vertex.h>
#include <core/cell/side.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/cell.h>
#include <core/generator/structured_decomposition.h>
#include <utils/math.h>
#include <core/face/face.h>
#include <core/generator/sector3d_generator.h>
#include <utils/geom.h>

#if DIM3

Sector3dGenerator::Sector3dGenerator(const Configuration &config)
        : StructuredGenerator(GeometryType::SECTOR) {

    m_global_r1 = config("mesh", "geometry", "R1").to_double();
    m_global_r2 = config("mesh", "geometry", "R2").to_double();

    m_alpha = config("mesh", "geometry", "Alpha").to_double();
    m_beta  = config("mesh", "geometry", "Beta").to_double();

    // В секторе количество ячеек по Y определяется алгоритмически
    m_global_nx = config("mesh", "cells_per_x").to_uint();
    m_global_ny = (size_t)floor(m_beta * m_global_nx / m_alpha) + 1;
    if (config.exist("mesh", "cells_per_y")) {
        m_global_ny = config("mesh", "cells_per_y").to_uint();
    }
    m_global_nz = static_cast<size_t>(math::int_log(m_global_r2 / m_global_r1, 1 + m_alpha / m_global_nx));

    m_type = GeometryType::SECTOR3D;
}

vector<vector<vector<Vertex::Ptr>>> Sector3dGenerator::create_vertices(Decomposition *decomp) const {
    double r1 = m_global_r1;
    double r2 = m_global_r2;

    double alpha_1 = - 0.5 * m_alpha;
    double alpha_2 = + 0.5 * m_alpha;
    double beta_1  = - 0.5 * m_beta;
    double beta_2  = + 0.5 * m_beta;

    double d_alpha = (alpha_2 - alpha_1) / m_global_nx;
    double d_beta  = (beta_2 - beta_1) / m_global_ny;

    auto vertices = vector<vector<vector<Vertex::Ptr>>>(m_global_nx + 1);
    for (size_t i = 0; i <= m_global_nx; ++i) {
        double phi1 = alpha_1 + d_alpha * i;
        double tg_phi1 = tan(phi1);

        vertices[i] = vector<vector<Vertex::Ptr>>(m_global_ny + 1);
        for (size_t j = 0; j <= m_global_ny; ++j) {
            double phi2 = beta_1 + d_beta * j;
            double tg_phi2 = tan(phi2);

            vertices[i][j] = vector<Vertex::Ptr>(m_global_nz + 1);
            for (size_t k = 0; k <= m_global_nz; ++k) {
                double r;

                if (k == m_global_nz) {
                    r = r2;
                } else if (k == 0) {
                    r = r1;
                } else {
                    r = r1 * pow(1 + d_alpha, k);
                }


                double z = r / sqrt(1.0 + tg_phi1 * tg_phi1 + tg_phi2 * tg_phi2);

                double x = z * tg_phi1;
                double y = z * tg_phi2;

                vertices[i][j][k] = Vertex::create(x, y, z);
            }
        }
    }
    return vertices;
}    

Vector3d loc2cart(double r, double tg1, double tg2) {
    Vector3d res = {0.0, 0.0, 0.0};
    res[2] = r / sqrt(1 + tg1 * tg1 + tg2 * tg2);
    res[0] = r * tg1;
    res[1] = r * tg2;
    return res;
};

Cell::Ptr Sector3dGenerator::create_base_cell(size_t u) const {
    auto p = u2xyz(u);
    size_t x = p[0];
    size_t y = p[1];
    size_t z = p[2];

    double da = m_alpha / m_global_nx;
    double db = m_beta  / m_global_ny;
    
    double a1 = -0.5 * m_alpha + da * x;
    double a2 = -0.5 * m_alpha + da * (x + 1);    
    double b1 = -0.5 * m_beta + db * y;
    double b2 = -0.5 * m_beta + db * (y + 1);
    
    double r1 = m_global_r1 * pow(1 + da, z);
    double r2 = z < m_global_nz ? m_global_r1 * pow(1 + da, z + 1) : m_global_r2;
    
    double tg_a1 = tan(a1);
    double tg_a2 = tan(a2);
    double tg_b1 = tan(b1);
    double tg_b2 = tan(b2);

    array<Vertex::Ptr, 8> v = {
            Vertex::create(loc2cart(r1, tg_a1, tg_b1)),
            Vertex::create(loc2cart(r1, tg_a2, tg_b1)),
            Vertex::create(loc2cart(r1, tg_a2, tg_b2)),
            Vertex::create(loc2cart(r1, tg_a1, tg_b2)),
            Vertex::create(loc2cart(r2, tg_a1, tg_b1)),
            Vertex::create(loc2cart(r2, tg_a2, tg_b1)),
            Vertex::create(loc2cart(r2, tg_a2, tg_b2)),
            Vertex::create(loc2cart(r2, tg_a1, tg_b2))
    };

    auto left_face   = Face::create(v[0], v[4], v[5], v[1]);
    auto right_face  = Face::create(v[3], v[7], v[6], v[2]);
    auto bottom_face = Face::create(v[0], v[1], v[2], v[3]);
    auto top_face    = Face::create(v[4], v[5], v[6], v[7]);
    auto back_face   = Face::create(v[0], v[3], v[7], v[4]);
    auto front_face  = Face::create(v[1], v[2], v[6], v[7]);

    auto cell = Cell::create(left_face, bottom_face, right_face, top_face, front_face, back_face);
    cell->set_z(xyz2u(x, y, z));
    cell->set_id(cell->z());

    return cell;
}

inline double sqr(double x) {
    return x * x;
}

Face::Ptr Sector3dGenerator::next_face_by_side(Cell::Ref cell, Side side) const {
    if (cell->faces(side)->is_complex()) {
        throw runtime_error("Sector3dGenerator::next_face_by_side(...) error: Complex face");
    }

    auto face = cell->faces(side)->single();

    Vector3d f_center = face->center();
    Vector3d opp_center = cell->faces(opposite_side(side))->center();

    std::function<Vertex::Ptr(Vertex::Ref)> transform;

    if (side == Side::LEFT || side == Side::RIGHT) {
        double tg_a0 = opp_center[0] / opp_center[2];
        double tg_a1 = f_center[0] / f_center[2];

        double tg_a2 = (2 * tg_a1 - tg_a0 * (1 - sqr(tg_a1))) / (1 - sqr(tg_a1) + 2 * tg_a0 * tg_a1);

        transform = [tg_a1, tg_a2](Vertex::Ref vertex) -> Vertex::Ptr {
            Vector3d v = vertex->v();

            double tg_b = v[1] / v[2];

            double x, y, z;
            z = sqrt(v.squaredNorm()/ (1.0 + sqr(tg_b) + sqr(tg_a2)));
            x = z * tg_a2;
            y = z * tg_b;
            return Vertex::create(x, y, z);
        };
    }
    else if (side == Side::BOTTOM || side == Side::TOP) {
        double tg_b0 = opp_center[1] / opp_center[2];
        double tg_b1 = f_center[1] / f_center[2];

        double tg_b2 = (2 * tg_b1 - tg_b0 * (1 - sqr(tg_b1))) / (1 - sqr(tg_b1) + 2 * tg_b0 * tg_b1);

        transform = [tg_b1, tg_b2](Vertex::Ref vertex) -> Vertex::Ptr {
                    Vector3d v = vertex->v();

                    double tg_a = v[0] / v[2];

                    double x, y, z;
                    z = sqrt(v.squaredNorm()/ (1.0 + sqr(tg_a) + sqr(tg_b2)));
                    x = z * tg_a;
                    y = z * tg_b2;
                    return Vertex::create(x, y, z);
                };
    }
    else {
        double r0 = opp_center.norm();
        double r1 = f_center.norm();
        double r2 = 2 * r1 - r0;

        transform = [r1, r2](Vertex::Ref vertex) -> Vertex::Ptr {
            double x = vertex->x() * r2 / r1;
            double y = vertex->y() * r2 / r1;
            double z = vertex->z() * r2 / r1;
            return Vertex::create(x, y, z);
        };
    }

    return Face::create(
            transform(face->vertex(0)),
            transform(face->vertex(1)),
            transform(face->vertex(2)),
            transform(face->vertex(3))
    );
}

void Sector3dGenerator::print_info(const std::string& tab) const {
    std::cout << tab << "Inner radius: " << std::hash<double>()(m_global_r1) << "\n";
    std::cout << tab << "Outer radius: " << std::hash<double>()(m_global_r2) << "\n";
    std::cout << tab << "Angle1:       " << std::hash<double>()(m_alpha) << "\n";
    std::cout << tab << "Angle2:       " << std::hash<double>()(m_beta) << "\n";
}

double Sector3dGenerator::r_in() const {
    return m_global_r1;
}

double Sector3dGenerator::r_out() const {
    return m_global_r2;
}

double Sector3dGenerator::angle1() const {
    return m_alpha;
}

double Sector3dGenerator::angle2() const {
    return m_beta;
}

#endif