//
// Created by 159-mrv on 3/28/18.
//

#include <control/mpi_wrapper.h>
#include <control/configuration.h>

#include <core/cell/side.h>
#include <core/generator/structured_decomposition.h>

StructuredDecomposition::StructuredDecomposition(const Configuration &config, array<size_t, 2> sizes) {
    read_config(config);
    processes_decomposition(sizes[0], sizes[1]);
    mesh_decomposition(sizes[0], sizes[1]);
}

void StructuredDecomposition::read_config(const Configuration &config) {
    if (config.exist("mesh", "decomposition")) {
        auto dec = config("mesh", "decomposition").to_string();

        if (dec == "X") {
            m_type = Type::X;
        } else if (dec == "Y") {
            m_type = Type::Y;
        } else if (dec == "XY") {
            m_type = Type::XY;
        } else {
            throw runtime_error("Unknown decomposition type, available values: 'X', 'Y', 'XY'");
        }
    } else {
        throw runtime_error("Configuration file doesn't contain mesh.decomposition parameter!");
    }

    m_pr_nx = m_pr_ny = -1;
    if (config.exist("mesh", "proc_per_x")) {
        m_pr_nx = config("mesh", "proc_per_x").to_int();
    }
    if (config.exist("mesh", "proc_per_y")) {
        m_pr_ny = config("mesh", "proc_per_y").to_int();
    }
}

void StructuredDecomposition::processes_decomposition(size_t nx, size_t ny) {
    if (m_type != Type::XY) {
        processes_decomposition_1d(static_cast<int>(nx), static_cast<int>(ny));
    } else {
        processes_decomposition_2d(static_cast<int>(nx), static_cast<int>(ny));
    }

    m_pr_x = mpi::rank() % m_pr_nx;
    m_pr_y = mpi::rank() / m_pr_nx;
}

void StructuredDecomposition::processes_decomposition_1d(int nx, int ny) {
    if (m_pr_nx >= 0 || m_pr_ny >= 0) {
        if (mpi::is_master()) {
            std::cerr << "  Warning: Configuration file includes parameter 'mesh.proc_per_x' or/and 'mesh.proc_per_y'\n";
            std::cerr << "           which has no sense for one dimensional decomposition.\n";
            std::cerr << "           The values will be ignored.\n\n";
        }
    }
    if (m_type == Type::X) {
        if (nx < mpi::size()) {
            throw runtime_error("Number of cells per X is less than number of processes.");
        }
        m_pr_nx = mpi::size();
        m_pr_ny = 1;
    } else {
        if (ny < mpi::size()) {
            throw runtime_error("Number of cells per Y is less than number of processes.");
        }
        m_pr_nx = 1;
        m_pr_ny = mpi::size();
    }
}

void StructuredDecomposition::processes_decomposition_2d(int nx, int ny) {
    if (mpi::size() == 1) {
        m_pr_nx = 1;
        m_pr_ny = 1;
        return;
    }

    if (m_pr_nx >= 0 && m_pr_ny >= 0) {
        if (m_pr_nx * m_pr_ny == mpi::size()) {
            return;
        }
        else {
            if (mpi::is_master()) {
                std::cerr << "  Warning: Configuration file includes parameters 'mesh.proc_per_x' and 'mesh.proc_per_y',\n";
                std::cerr << "           but proc_per_x * proc_per_y != mpi::size().\n";
                std::cerr << "           The values will be ignored.\n\n";
            }
        }
    }

    if (m_pr_nx >= 0) {
        if (mpi::size() % m_pr_nx == 0) {
            m_pr_ny = mpi::size() / m_pr_nx;
            return;
        }
        else {
            if (mpi::is_master()) {
                std::cerr << "  Warning: Configuration file includes parameters 'mesh.proc_per_x',\n";
                std::cerr << "           but mpi::size() is not a multiple of the value.\n";
                std::cerr << "           The value will be ignored.\n\n";
            }
        }
    }

    if (m_pr_ny >= 0) {
        if (mpi::size() % m_pr_ny == 0) {
            m_pr_nx = mpi::size() / m_pr_ny;
            return;
        }
        else {
            if (mpi::is_master()) {
                std::cerr << "  Warning: Configuration file includes parameters 'mesh.proc_per_y',\n";
                std::cerr << "           but mpi::size() is not a multiple of the value.\n";
                std::cerr << "           The value will be ignored.\n\n";
            }
        }
    }


    // Находим делители
    vector<int> l1, l2;

    for (int n = 1; n <= mpi::size(); ++n) {
        if (n > nx || mpi::size() / n > ny) {
            continue;
        }
        if (mpi::size() % n == 0) {
            l1.push_back(n);
            l2.push_back(mpi::size() / n);
        }
    }

#if 0
    // В данном блоке минимизируется функционал
    // F(pr_nx, pr_ny) = |pr_nx - pr_ny|
    // Таким образом, выбирается близкое число процессов по осям Ox и Oy,
    // то есть pr_nx ~ pr_ny.
    // Так, если число процессов равно 12, то получим разбиение 3 x 4,
    // независимо от размеров области

    int min_diff = std::numeric_limits<int>::max();;
    for (size_t i = 0; i < l1.size(); ++i) {
        auto diff = abs(l1[i] - l2[i]);

        if (min_diff > diff) {
            m_pr_nx = l1[i];
            m_pr_ny = l2[i];

            min_diff = diff;
        }
    }
#else
    // В данном блоке минимизируется функционал
    // F(pr_nx, pr_ny) = |ny * pr_nx - nx pr_ny|
    // Области разбиения выбираются так, чтобы они были близки к квадратам,
    // то есть nx / pr_nx ~ ny / pr_ny.
    // Так, если число процессов равно 12, а область вытянута вдоль оси Ox,
    // то можно ожидать разбиение типа 4 x 3, 6 x 2, 12 x 1.

    int min_diff = std::numeric_limits<int>::max();
    for (size_t i = 0; i < l1.size(); ++i) {
        auto diff = abs(ny * l1[i] - nx * l2[i]);

        if (min_diff > diff) {
            m_pr_nx = l1[i];
            m_pr_ny = l2[i];

            min_diff = diff;
        }
    }
#endif
}

void StructuredDecomposition::mesh_decomposition(size_t nx, size_t ny) {
    // Размер небольшого блока (большой на единицу больше)
    size_t x_block_size = nx / m_pr_nx;
    size_t y_block_size = ny / m_pr_ny;

    // Количество больших блоков по осям X и Y
    int x_big_num = static_cast<int>(nx) % m_pr_nx;
    int y_big_num = static_cast<int>(ny) % m_pr_ny;

    // Является ли текущий блок большим
    bool x_is_big = m_pr_x < x_big_num;
    bool y_is_big = m_pr_y < y_big_num;

    m_nx = x_block_size + static_cast<size_t>(x_is_big);
    m_ny = y_block_size + static_cast<size_t>(y_is_big);

    m_x = x_block_size * m_pr_x + (x_is_big ? m_pr_x : x_big_num);
    m_y = y_block_size * m_pr_y + (y_is_big ? m_pr_y : y_big_num);
}

int StructuredDecomposition::cell_rank(size_t x, size_t y) const noexcept {
    int _pr_x = m_pr_x;
    int _pr_y = m_pr_y;

    if (x < m_x) {
        --_pr_x;
    }
    if (x >= m_x + m_nx) {
        ++_pr_x;
    }
    if (y < m_y) {
        --_pr_y;
    }
    if (y >= m_y + m_ny) {
        ++_pr_y;
    }

    return _pr_y * m_pr_nx + _pr_x;
}