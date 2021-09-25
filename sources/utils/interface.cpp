//
// Created by 159-mrv on 7/17/19.
//

#include <control/configuration.h>

#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/cell/cell.h>

#include <io/vtk/vtk_writer.h>

#include <utils/memory/node_list.h>
#include <utils/geom.h>
#include <utils/math.h>

#include <utils/interface.h>
#include <control/mpi_wrapper.h>


Interface::Interface() :
    m_size_x(0), m_size_y(0), m_r_min(0.0), m_r_max(0.0), m_h2(0.0), m_length(0.0) {

}

Interface::~Interface() = default;

Interface::Interface(const vector<vector<Vector3d>>& points) :
    m_size_x(points.size()),
    m_size_y(points.front().size()),
    m_points(points) {
    Vector3d inf = {math::inf(), math::inf(), math::inf()};
    m_offsets = vector<vector<Vector3d>>(m_size_x, vector<Vector3d>(m_size_y, inf));
    m_dists2 = vector<vector<double>>(m_size_x, vector<double>(m_size_y, math::inf()));

    calc_params();
}

void Interface::write(string filename) const {
    NodeList::Ptr cells = NodeList::create(0);

    double xi = 1.0e-3 * m_length;

#if DIM2
    for (size_t i = 1; i < m_size_x; ++i) {
        Vector3d normal = {m_points[i][0][1] - m_points[i - 1][0][1],
                           m_points[i - 1][0][0] - m_points[i][0][0],
                           0.0};
        normal.normalize();
        normal = xi * normal;

        Vertex::Ptr v1 = Vertex::create(m_points[i - 1][0] - normal);
        Vertex::Ptr v2 = Vertex::create(m_points[i - 1][0] + normal);
        Vertex::Ptr v3 = Vertex::create(m_points[i][0] + normal);
        Vertex::Ptr v4 = Vertex::create(m_points[i][0] - normal);

        Face::Ptr left = Face::create(v1, v2);
        Face::Ptr right = Face::create(v3, v4);
        Face::Ptr top = Face::create(v2, v3);
        Face::Ptr bottom = Face::create(v1, v4);

        Cell::Ptr cell = Cell::create(left, bottom, right, top);
        cells->push_back(cell);
    }
#else
    for (size_t i = 1; i < m_size_x; ++i) {
        for (size_t j = 1; j < m_size_y; ++j) {
            array<Vector3d, 4> bv = {
                    m_points[i][j], m_points[i-1][j], m_points[i-1][j-1], m_points[i][j-1]
            };

            Vector3d normal = geom::normal(bv);
            normal = xi * normal;

            array<Vertex::Ptr, 8> v = {
                    Vertex::create(bv[0] - normal),
                    Vertex::create(bv[1] - normal),
                    Vertex::create(bv[2] - normal),
                    Vertex::create(bv[3] - normal),
                    Vertex::create(bv[0] + normal),
                    Vertex::create(bv[1] + normal),
                    Vertex::create(bv[2] + normal),
                    Vertex::create(bv[3] + normal)
            };

            Face::Ptr left   = Face::create(v[0], v[3], v[7], v[4]);
            Face::Ptr right  = Face::create(v[1], v[2], v[6], v[5]);
            Face::Ptr bottom = Face::create(v[0], v[1], v[5], v[4]);
            Face::Ptr top    = Face::create(v[7], v[6], v[2], v[3]);
            Face::Ptr back   = Face::create(v[0], v[1], v[2], v[3]);
            Face::Ptr front  = Face::create(v[4], v[5], v[6], v[7]);

            Cell::Ptr cell = Cell::create(left, bottom, right, top, front, back);
            cells->push_back(cell);
        }
    }
#endif

    VtkWriter::write(*cells, filename);
}

void Interface::new_line_vertex(const Vector3d& v) {
    for(size_t i = 0; i < m_size_x; ++i) {
        for(size_t j = 0; j < m_size_y; ++j) {
            double dist2 = (m_points[i][j] - v).squaredNorm();
            if (dist2 < m_dists2[i][j]) {
                m_offsets[i][j] = v - m_points[i][j];
                m_dists2[i][j] = dist2;
            }
        }
    }
}

void Interface::move(bool smooth) {
    // Крайние значения проецируются на границы сектора
    for (size_t i: {size_t(0), m_size_x - 1}) {
        for (size_t j = 0; j < m_size_y; ++j) {
            Vector3d off = m_offsets[i][j];
            Vector3d poi = m_points[i][j];

            double a = poi.squaredNorm();
            double b =  off.dot(poi);

            Vector3d keks = m_points[i][j] * b / a;
            m_offsets[i][j] = keks;
        }
    }
    if (m_size_y > 1) {
        for (size_t j: {size_t(0), m_size_y - 1}) {
            for (size_t i = 0; i < m_size_x; ++i) {
                m_offsets[i][j] = m_points[i][j] * m_offsets[i][j].dot(m_points[i][j]) / m_points[i][j].squaredNorm();
            }
        }
    }

    // Крайние значения проецируются на границы прямоугольника
    //m_offsets.front()[1] = 0.0;
    //m_offsets.back()[0] = 0.0;


    for(size_t i = 0; i < m_size_x; ++i) {
        for(size_t j = 0; j < m_size_y; ++j) {
            m_points[i][j] += m_offsets[i][j];
        }
    }

    if (smooth) {
        // Сглаживание
        for (size_t i = 1; i < m_size_x - 1; ++i) {
            for(size_t j = 1; j < m_size_y - 1; ++j) {
                m_offsets[i][j] = (m_points[i - 1][j] - 2 * m_points[i][j] + m_points[i + 1][j] +
                                   m_points[i][j - 1] - 2 * m_points[i][j] + m_points[i][j + 1]) / 8.0;
            }
        }

        for (size_t i = 1; i < m_size_x - 1; ++i) {
            for(size_t j = 1; j < m_size_y - 1; ++j) {
                m_points[i][j] += m_offsets[i][j];
            }
        }
    }

    // Зануляем
    for(size_t i = 0; i < m_size_x; ++i) {
        for(size_t j = 0; j < m_size_y; ++j) {
            m_offsets[i][j] = {0.0, 0.0, 0.0};
            m_dists2[i][j] = math::inf();
        }
    }

    calc_params();
}

double Interface::sqr_distance(const Vector3d &v) const {
    //double probe_dist = (v - m_points[m_size_x / 2][m_size_y / 2]).squaredNorm();
   // if (probe_dist > 0.09 * m_length * m_length) {
   //     return math::inf();
  //  }

    double min_dist2 = math::inf();
    if (m_size_y < 2) {
        for(size_t i = 1; i < m_size_x; ++i) {
            double dist2 = geom::sqr_distance(v, m_points[i][0], m_points[i - 1][0]);
            if (dist2 < min_dist2) {
                min_dist2 = dist2;
            }
        }
    }
    else {
        for (size_t i = 1; i < m_size_x; ++i) {
            for (size_t j = 1; j < m_size_y; ++j) {
                double dist2 = geom::sqr_distance(v, m_points[i][j], m_points[i - 1][j],
                                                  m_points[i - 1][j - 1], m_points[i][j - 1]);
                if (dist2 < min_dist2) {
                    min_dist2 = dist2;
                }
            }
        }
    }
    return min_dist2;
}

void Interface::calc_params() {
    //m_h2 = 0.02;
    //return;

    m_r_min = math::inf();
    m_r_max = 0.0;
    for (auto &row: m_points) {
        for (auto &v: row) {
            double r2 = v.squaredNorm();
            m_r_max = std::max(m_r_max, r2);
            m_r_min = std::min(m_r_min, r2);
        }
    }
    m_r_max = std::sqrt(m_r_max);
    m_r_min = std::sqrt(m_r_min);

    m_length = 0.0;
    for (size_t i = 1; i < m_size_x; ++i) {
        m_length += (m_points[i][0] - m_points[i - 1][0]).norm();
    }

    // for 1 rnp
    // m_h2 = 9.0e-4 * m_length * m_length;

    // for 1 shrh
    //m_h2 = 4.0 * m_length * m_length;

    // for 5 shrh
    m_h2 = 0.1 * m_length * m_length;
}

#ifdef ENABLE_MPI
void norm2_min(void* _in, void* _inout, int* len, MPI_Datatype *dptr) {
    Vector3d* in = (Vector3d*)_in;
    Vector3d* inout = (Vector3d*)_inout;

    for(int i = 0; i < *len; ++i) {
        double n1 = in[i].squaredNorm();
        double n2 = inout[i].squaredNorm();

        if (n1 < n2) {
            inout[i] = in[i];
        }
    }
}
#endif

void Interface::reduce(const vector<Interface>& bords) {
    size_t num = bords.size();

    vector<double> min_dist2(m_size_x * m_size_y);
    vector<Vector3d> min_off(m_size_x * m_size_y);

    // минимизация по тредам
    for(size_t i = 0; i < m_size_x; ++i) {
        for (size_t j = 0; j < m_size_y; ++j) {
            min_dist2[m_size_y * i + j] = bords.front().m_dists2[i][j];
            size_t min_idx = 0;
            for (size_t k = 1; k < num; ++k) {
                if (bords[k].m_dists2[i][j] < min_dist2[m_size_y * i + j]) {
                    min_dist2[m_size_y * i + j] = bords[k].m_dists2[i][j];
                    min_idx = k;
                }
            }
            min_off[m_size_y * i + j] = bords[min_idx].m_offsets[i][j];
        }
    }

#ifdef ENABLE_MPI
    MPI_Op MPI_NORM2_MIN;
    MPI_Datatype MPI_VECTOR3D;

    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_VECTOR3D);
    MPI_Type_commit(&MPI_VECTOR3D);
    MPI_Op_create(&norm2_min, true, &MPI_NORM2_MIN);

    vector<Vector3d> offsets = vector<Vector3d>(m_size_x * m_size_y);

    MPI_Allreduce(min_off.data(), offsets.data(), int(m_size_x*m_size_y), MPI_VECTOR3D, MPI_NORM2_MIN, MPI_COMM_WORLD);

    for(size_t i = 0; i < m_size_x; ++i) {
        for(size_t j = 0; j < m_size_y; ++j) {
            m_offsets[i][j] = offsets[m_size_y * i + j];
        }
    }
#else
    for(size_t i = 0; i < m_size_x; ++i) {
        for(size_t j = 0; j < m_size_y; ++j) {
            m_offsets[i][j] = min_off[m_size_y * i + j];
        }
    }
#endif
}
