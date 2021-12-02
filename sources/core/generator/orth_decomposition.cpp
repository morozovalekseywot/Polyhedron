//
// Created by 159-mrv on 2/15/19.
//

#include <control/configuration.h>
#include <core/generator/orth_decomposition.h>
#include <utils/geom.h>
#include <control/mpi_wrapper.h>
#include <core/generator/geometry_generator.h>
#include <core/generator/real_sector_generator.h>
#include <core/generator/sector_generator.h>
#include <core/generator/rectangle_generator.h>
#include <core/generator/cuboid_generator.h>
#include <utils/math.h>
#include <core/generator/sector3d_generator.h>
#include <random>


OrthDecomposition::OrthDecomposition(const Configuration &config, GeometryGenerator* geometry) {
    m_type = geometry->type();

    read_config(config);

    auto c = proc_coord(mpi::rank());
    m_x = c[0];
    m_y = c[1];
    m_z = c[2];

    init_barriers(geometry);

    skew_barriers();

    m_dx = vector<double>(size_t(m_px - 1), 0.0);
    m_dy = vector<vector<double>>(size_t(m_px));
    m_dz = vector<vector<vector<double>>>(size_t(m_px));
    for (int i = 0; i < m_px; ++i) {
        m_dy[i] = vector<double>(size_t(m_py[i] - 1), 0.0);
        m_dz[i] = vector<vector<double>>(size_t(m_py[i]));
        for (int j = 0; j < m_py[i]; ++j) {
            m_dz[i][j] = vector<double>(size_t(m_pz[i][j] - 1), 0.0);
        }
    }

    m_workload = vector<vector<vector<double>>>(size_t(m_px));
    for (int i = 0; i < m_px; ++i) {
        m_workload[i] = vector<vector<double>>(size_t(m_py[i]));
        for (int j = 0; j < m_py[i]; ++j) {
            m_workload[i][j] = vector<double>(size_t(m_pz[i][j]), 0.0);
        }
    }
}

void OrthDecomposition::read_config(const Configuration& config) {
    if (mpi::size() == 1) {
        m_dir = {Direction::Ox, Direction::Oy, Direction::Oz};

        m_px = 1;
        m_py = {1};
        m_pz = {{1}};

        return;
    }

    static_vector<Direction, 3> dirs;
    if (config.exist("mesh", "decomposition")) {
        string dir = config("mesh", "decomposition").to_string();
        for (size_t i = 0; i < min(dir.size(), size_t(3)); ++i) {
            switch (dir[i]) {
                case 'X':
                    dirs.push_back(Direction::Ox);
                    break;
                case 'Y':
                    dirs.push_back(Direction::Oy);
                    break;
                case 'Z':
                    dirs.push_back(Direction::Oz);
                    break;
                default:
                    throw runtime_error("Set correct 'decomposition' in config");
            }
        }
    } else {
        throw runtime_error("Set parameter 'decomposition' in config");
    }

    for (auto d1: {Direction::Ox, Direction::Oy, Direction::Oz}) {
        bool contain = false;
        for (auto d2: dirs) {
            if (d1 == d2) {
                contain = true;
                break;
            }
        }
        if (!contain) {
            dirs.push_back(d1);
        }
    }
    m_dir = {dirs[0], dirs[1], dirs[2]};

    static_vector<uint, 3> proc;
    for (auto dir: m_dir) {
        string proc_per = "proc_per_";
        switch (dir) {
            case Direction::Ox:
                proc_per += "x";
                break;
            case Direction::Oy:
                proc_per += "y";
                break;
            case Direction::Oz:
                proc_per += "z";
                break;
        }

        uint pr = 1;
        if (config.exist("mesh", proc_per)) {
            pr = config("mesh", proc_per).to_uint();
        }
        proc.push_back(pr);
    }

    if (proc[0] * proc[1] * proc[2] != mpi::usize()) {
        throw runtime_error("Decomposition sizes are inconsistent");
    }

    m_px = proc[0];
    m_py = vector<int>(proc[0], proc[1]);
    m_pz = vector<vector<int>>(proc[0], vector<int>(proc[1], proc[2]));


#if DIM2    // Ограничения двумерной геометрии
    if (m_dir[0] == Direction::Oz || m_dir[1] == Direction::Oz) {
        throw runtime_error(
                "OrthDecomposition error: 'Z' decomposition is not avaliable for 2D domain");
    }

    if (m_pz[0][0] > 1) {
        throw runtime_error("OrthDecomposition error: 'proc_per_z' > 1 for 2D calculation domain");
    }
#endif
}

void OrthDecomposition::init_barriers(GeometryGenerator* geometry) {
    switch (m_type) {
#if DIM2
        case GeometryType::SECTOR:
            sector_barriers_init((SectorGenerator*)geometry);
            break;

        case GeometryType::REAL_SECTOR:
            real_sector_barriers_init((RealSectorGenerator*)geometry);
            break;

        case GeometryType::RECTANGLE:
            rectangle_barriers_init((RectangleGenerator*)geometry);
            break;
#else
        case GeometryType::SECTOR3D:
            sector3d_barriers_init((Sector3dGenerator*)geometry);
            break;

        case GeometryType::CUBOID:
            cuboid_barriers_init((CuboidGenerator*)geometry);
            break;
#endif
        default:
            throw runtime_error("Strange geometry for OrthDecomposition");

    }
}

void OrthDecomposition::sector_barriers_init(SectorGenerator* sector) {
    double r_0 = sector->r_in();
    double r_in = sector->r_in();
    double r_out = sector->r_out();
    double alpha = sector->angle();

    common_sector_barriers_init(r_0, r_in, r_out, alpha);
}

void OrthDecomposition::real_sector_barriers_init(RealSectorGenerator* sector) {
    double r_0 = 0.0;
    double r_in = sector->r_min();
    double r_out = sector->r_out();
    double alpha = sector->angle();

    common_sector_barriers_init(r_0, r_in, r_out, alpha);
}

void OrthDecomposition::common_sector_barriers_init(double r0, double r1, double r2, double alpha) {
    m_to_local = [this](const Vector3d &v) -> Vector3d {
        Vector3d v1 = {
                atan2(v[0], v[1]),
                sqrt(v[0] * v[0] + v[1] * v[1]),
                v[2] };
        return {
                v1[static_cast<int>(m_dir[0])],
                v1[static_cast<int>(m_dir[1])],
                v1[static_cast<int>(m_dir[2])]
        };
    };

    auto get_coord = [r0, r1, r2, alpha](Direction dir, double xi) -> double {
        if (dir == Direction::Ox) {
            return (xi - 0.5) * alpha;

        } else if (dir == Direction::Oy) {
            //return xi > 0.0 ? r1 * pow(r2 / r1, xi) : r0;

            double l0 = r0 / r1;
            double l2 = r2 / r1;

            double k = M_PI / 4.0;

            double K = sqrt((1 - xi) * l0 * l0 + xi * (1 + (2.0 / k) * log(l2)));
            if (K > 1.0) {
                K = exp(0.5 * k * (xi - 1) * (1 - l0 * l0)) * pow(l2, xi);
            }
            return r1 * K;

        } else {
            return (xi - 0.5) * math::inf();
        }
    };

    set_barriers(get_coord);
}

void OrthDecomposition::rectangle_barriers_init(RectangleGenerator* rect) {
    m_to_local = [this](const Vector3d &v) -> Vector3d {
        return {
                v[static_cast<int>(m_dir[0])],
                v[static_cast<int>(m_dir[1])],
                v[static_cast<int>(m_dir[2])]
        };
    };

    double x_min = rect->x_min();
    double x_len = rect->x_len();
    double y_min = rect->y_min();
    double y_len = rect->y_len();

    auto get_coord = [x_min, x_len, y_min, y_len](Direction dir, double xi) -> double {
        if (dir == Direction::Ox) {
            return x_min + xi * x_len;

        } else if (dir == Direction::Oy) {
            return y_min + xi * y_len;

        } else {
            return (xi - 0.5) * math::inf();
        }
    };

    set_barriers(get_coord);
}

#if DIM3

void OrthDecomposition::sector3d_barriers_init(Sector3dGenerator* sector) {
    m_to_local = [this](const Vector3d &v) -> Vector3d {
        Vector3d v1 = {
                atan2(v[0], v[2]),
                atan2(v[1], v[2]),
                v.norm()
        };
        return {
                v1[static_cast<int>(m_dir[0])],
                v1[static_cast<int>(m_dir[1])],
                v1[static_cast<int>(m_dir[2])]
        };
    };

    double r1 = sector->r_in();
    double r2 = sector->r_out();
    double alpha = sector->angle1();
    double beta = sector->angle2();

    auto get_coord = [r1, r2, alpha, beta](Direction dir, double xi) -> double {
        if (dir == Direction::Ox) {
            return (xi - 0.5) * alpha;

        } else if (dir == Direction::Oy) {
            return (xi - 0.5) * beta;

        } else {
            return r1 * pow(r2 / r1, xi);
        }
    };

    set_barriers(get_coord);
}

void OrthDecomposition::cuboid_barriers_init(CuboidGenerator* cuboid) {
    m_to_local = [this](const Vector3d &v) -> Vector3d {
        return {
                v[static_cast<int>(m_dir[0])],
                v[static_cast<int>(m_dir[1])],
                v[static_cast<int>(m_dir[2])]
        };
    };

    double x_min = cuboid->x_min();
    double x_len = cuboid->x_len();
    double y_min = cuboid->y_min();
    double y_len = cuboid->y_len();
    double z_min = cuboid->z_min();
    double z_len = cuboid->z_len();

    auto get_coord = [x_min, x_len, y_min, y_len, z_min, z_len](Direction dir, double xi) -> double {
        if (dir == Direction::Ox) {
            return x_min + xi * x_len;

        } else if (dir == Direction::Oy) {
            return y_min + xi * y_len;

        } else {
            return z_min + xi * z_len;
        }
    };

    set_barriers(get_coord);
}

#endif

void OrthDecomposition::set_barriers(const std::function<double(Direction, double)> &get_coord) {
    m_barriers_x = vector<double>(size_t(m_px) + 1);
    m_barriers_y = vector<vector<double>>(size_t(m_px));
    m_barriers_z = vector<vector<vector<double>>>(size_t(m_px));

    for (int i = 0; i < m_px; ++i) {
        m_barriers_x[i] = get_coord(m_dir[0], double(i) / m_px);

        m_barriers_y[i] = vector<double>(size_t(m_py[i]) + 1);
        m_barriers_z[i] = vector<vector<double>>(size_t(m_py[i]));

        for (int j = 0; j < m_py[i]; ++j) {
            m_barriers_y[i][j] = get_coord(m_dir[1], double(j) / m_py[i]);

            m_barriers_z[i][j] = vector<double>(size_t(m_pz[i][j]) + 1);

            for(int k = 0; k <= m_pz[i][j]; ++k) {
                m_barriers_z[i][j][k] = get_coord(m_dir[2], double(k) / m_pz[i][j]);
            }
        }
        m_barriers_y[i][m_py[i]] = get_coord(m_dir[1], 1.0);
    }
    m_barriers_x[m_px] = get_coord(m_dir[0], 1.0);
}

void OrthDecomposition::skew_barriers() {
    // В данной функции "ломаем" идеальное разбиение, поскольку оно может
    // приводить к ошибкам, связанным с симметрией.  Ошибка возникала при
    // точном совпадении центра ячеек с положением барьера.

    double len_x = m_barriers_x.back() - m_barriers_x[0];
    double len_y = m_barriers_y[0].back() - m_barriers_y[0][0];
    double len_z = m_barriers_z[0][0].back() - m_barriers_z[0][0][0];
    double eps = 1.0e-8;

    std::mt19937_64 gen(0);

    std::uniform_real_distribution<double> rand_x(-eps * len_x, eps * len_x);
    std::uniform_real_distribution<double> rand_y(-eps * len_y, eps * len_y);
    std::uniform_real_distribution<double> rand_z(-eps * len_z, eps * len_z);

    for (int i = 1; i < m_px; ++i) {
        m_barriers_x[i] = m_barriers_x[i] + rand_x(gen);
    }

    for (int i = 0; i < m_px; ++i) {
        for (int j = 1; j < m_py[i]; ++j) {
            m_barriers_y[i][j] += rand_y(gen);
        }
    }

    for (int i = 0; i < m_px; ++i) {
        for (int j = 0; j < m_py[i]; ++j) {
            for (int k = 1; k < m_pz[i][j]; ++k) {
                m_barriers_z[i][j][k] += rand_z(gen);
            }
        }
    }
}

array<int, 3> OrthDecomposition::proc_coord(int rank) const {
    return  {
            rank / (m_pz[0][0] * m_py[0]),
            (rank / m_pz[0][0]) % m_py[0],
            rank % m_pz[0][0]
    };
}

int OrthDecomposition::rank(const Vector3d& vertex) const {
    Vector3d v = m_to_local(vertex);

    //double _r   = v.norm();
    //double _phi = atan2(v[1], v[0]) - ROTATION_ANGLE;
    //v = {_r * cos(_phi), _r * sin(_phi), 0.0};

    // Определяем угол
    int I = -1;
    for (int i = 0; i < m_px; ++i) {
        if (m_barriers_x[i] <= v[0] && v[0] <= m_barriers_x[i + 1]) {
            I = i;
            break;
        }
    }
    if (I < 0) {
        return -1;
    }

    // Определяем радиус
    int J = -1;
    for (int j = 0; j < m_py[I]; ++j) {
        if (m_barriers_y[I][j] <= v[1] && v[1] <= m_barriers_y[I][j + 1]) {
            J = j;
            break;
        }
    }
    if (J < 0) {
        return -1;
    }

    // Третья координата
    int K = -1;
    for (int k = 0; k < m_pz[I][J]; ++k) {
        if (m_barriers_z[I][J][k] <= v[2] && v[2] <= m_barriers_z[I][J][k + 1]) {
            K = k;
            break;
        }
    }

    if (K < 0) {
        return -1;
    }

    return m_pz[0][0] * (m_py[0] * I + J) + K;
}

void OrthDecomposition::set_workload(double workload) {
    m_workload[m_x][m_y][m_z] = workload;
}

void OrthDecomposition::sync_data() {
    auto ws = mpi::all_gather(m_workload[m_x][m_y][m_z]);

    for(int r = 0; r < mpi::size(); ++r) {
        auto c = proc_coord(r);
        int x = c[0];
        int y = c[1];
        int z = c[2];
        m_workload[x][y][z] = ws[r];
    }
}

void OrthDecomposition::update() {
    // Первый этап
    for (int i = 0; i < m_px - 1; ++i) {
        double Li(0.0), Lj(0.0);

        double weight = 0.0;
        for (int j = 0; j < m_py[i]; ++j) {
            for (int k = 0; k < m_pz[i][j]; ++k) {
                Li += m_workload[i][j][k];
                weight += 1.0;
            }
        }
        Li /= weight;

        weight = 0.0;
        for (int j = 0; j < m_py[i + 1]; ++j) {
            for (int k = 0; k < m_pz[i + 1][j]; ++k) {
                Lj += m_workload[i + 1][j][k];
                weight += 1.0;
            }
        }
        Lj /= weight;

        double diff = (Lj - Li) / (Li + Lj);

        double Dij = 0.1 * min(fabs(m_barriers_x[i + 2] - m_barriers_x[i + 1]),
                               fabs(m_barriers_x[i + 1] - m_barriers_x[i]));
        m_dx[i] = Dij * diff;
    }
    for (int i = 1; i < m_px; ++i) {
        m_barriers_x[i] += m_dx[i - 1];
    }

    // Второй этап
    for (int i = 0; i < m_px; ++i) {
        for (int j = 0; j < m_py[i] - 1; ++j) {
            double Li(0.0), Lj(0.0);
            double weight = 0.0;
            for (int k = 0; k < m_pz[i][j]; ++k) {
                Li += m_workload[i][j][k];
                weight += 1.0;
            }
            Li /= weight;

            weight = 0.0;
            for (int k = 0; k < m_pz[i][j + 1]; ++k) {
                Lj += m_workload[i][j + 1][k];
                weight += 1.0;
            }
            Lj /= weight;

            double diff = (Lj - Li) / (Li + Lj);

            double Dij = 0.1 * min(fabs(m_barriers_y[i][j + 2] - m_barriers_y[i][j + 1]),
                                   fabs(m_barriers_y[i][j + 1] - m_barriers_y[i][j]));
            m_dy[i][j] = Dij * diff;
        }

        for (int j = 1; j < m_py[i]; ++j) {
            m_barriers_y[i][j] += m_dy[i][j - 1];
        }
    }

    // Третий этап
    for (int i = 0; i < m_px; ++i) {
        for (int j = 0; j < m_py[i]; ++j) {
            for (int k = 0; k < m_pz[i][j] - 1; ++k) {
                double Li = m_workload[i][j][k];
                double Lj = m_workload[i][j][k + 1];

                double diff = (Lj - Li) / (Li + Lj);

                double Dij = 0.1 * min(m_barriers_z[i][j][k + 2] - m_barriers_z[i][j][k + 1],
                                       m_barriers_z[i][j][k + 1] - m_barriers_z[i][j][k]);

                m_dz[i][j][k] = Dij * diff;
            }

            for (int k = 1; k < m_pz[i][j]; ++k) {
                m_barriers_z[i][j][k] += m_dz[i][j][k - 1];
            }
        }
    }
}
