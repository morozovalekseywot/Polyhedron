//
// Created by 159-mrv on 5/7/18.
//

#include <utils/partitioning/base_splitter_2d.h>
#include <control/mpi_wrapper.h>

BaseSplitter2D::BaseSplitter2D(const SplitterIPoint& sizes)
    : m_nx(sizes[0]), m_ny(sizes[1]), m_size(0) {
    m_workload = vector<vector<double>>(m_nx);
    m_rank = vector<vector<int>>(m_nx);
    for(size_t i = 0; i < m_nx; ++i) {
        m_workload[i] = vector<double>(m_ny, 0.0);
        m_rank[i] = vector<int>(m_ny, 33);
    }
}

void BaseSplitter2D::add_workload(const SplitterIPoint& point, double weight) {
    m_workload[point[0]][point[1]] += weight;
}

void BaseSplitter2D::add_workload(const SplitterFPoint &point, double weight) {
    m_workload[static_cast<size_t>(point[0])][static_cast<size_t>(point[1])] += weight;
}

void BaseSplitter2D::random_workload() {
    for(size_t i = 0; i < m_nx; ++i) {
        for(size_t j = 0; j < m_ny; ++j) {
            double x = static_cast<double>(i) / (m_nx - 1) - 0.75;
            double y = static_cast<double>(j) / (m_ny - 1) - 0.25;

            m_workload[i][j] = exp(-10 * (x * x + y * y));
        }
    }
}

void BaseSplitter2D::sync_data() {
    for(size_t i = 0; i < m_nx; ++i) {
        m_workload[i] = mpi::sum(m_workload[i]);
    }
}

void BaseSplitter2D::split(int size) {
    m_size = static_cast<uint>(size);
    Block b(0, 0, m_nx, m_ny);
    RankRange r(0, m_size);
    subsplit(b, r);
}

void BaseSplitter2D::subsplit(const Block &b, const RankRange &r, bool axis) {
    int n_rank = r.size();

    if (n_rank == 1) {
        for (size_t i = b.x1; i < b.x2; ++i) {
            for (size_t j = b.y1; j < b.y2; ++j) {
                m_rank[i][j] = r.begin;
            }
        }
        return;
    }

    // такс, мы здеся
    // определим в какую сторону разбивать
    auto nx = b.x_size();
    auto ny = b.y_size();

    bool split_x = nx > ny || (axis && nx == ny);

    vector<double> wls(split_x ? nx : ny, 0.0);
    vector<double> wl_sum(split_x ? nx : ny, 0.0);

    if (split_x) {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                wls[i] += m_workload[b.x1 + i][b.y1 + j];
            }
        }
    } else {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                wls[j] += m_workload[b.x1 + i][b.y1 + j];
            }
        }
    }


    wl_sum[0] = wls[0];
    for (size_t i = 1; i < wl_sum.size(); ++i) {
        wl_sum[i] = wl_sum[i - 1] + wls[i];
    }

    int min_range = n_rank / 2;
    int max_range = n_rank - min_range;


    double min_diff = 1.0 / 0.0;
    size_t delemeter = wls.size() * min_range / max_range;
    int left_rank_range = min_range;

    for (size_t i = 1; i < wls.size(); ++i) {
        double ratio = wl_sum[i - 1] / (wl_sum.back() - wl_sum[i - 1]);

        double diff1 = std::fabs(ratio - static_cast<double>(min_range) / max_range);
        double diff2 = std::fabs(ratio - static_cast<double>(max_range) / min_range);

        if (diff1 < min_diff) {
            min_diff = diff1;
            delemeter = i;
            left_rank_range = min_range;
        }

        if (diff2 < min_diff) {
            min_diff = diff2;
            delemeter = i;
            left_rank_range = max_range;
        }
    }

    Block left_block = split_x ? Block(b.x1, b.y1, b.x1 + delemeter, b.y2) :
                       Block(b.x1, b.y1, b.x2, b.y1 + delemeter);
    RankRange left_range(r.begin, r.begin + left_rank_range);
    subsplit(left_block, left_range, !axis);

    Block right_block = split_x ? Block(b.x1 + delemeter, b.y1, b.x2, b.y2) :
                        Block(b.x1, b.y1 + delemeter, b.x2, b.y2);
    RankRange right_range(r.begin + left_rank_range, r.end);
    subsplit(right_block, right_range, !axis);
}

void BaseSplitter2D::print_info(const string &tabs) const {
    vector<double> wl(m_size, 0.0);
    for (size_t i = 0; i < m_nx; ++i) {
        for (size_t j = 0; j < m_ny; ++j) {
            wl[m_rank[i][j]] += m_workload[i][j];
        }
    }

    print_workload_info(tabs, wl);
}

int BaseSplitter2D::rank(const SplitterIPoint& point) const {
    return m_rank[point[0]][point[1]];
}

int BaseSplitter2D::rank(const SplitterFPoint &point) const {
    return m_rank[static_cast<size_t>(point[0])][static_cast<size_t>(point[1])];
}

void BaseSplitter2D::print_ranks() const {
    for(size_t j = 0; j < m_ny; ++j) {
        for(size_t i = 0; i < m_nx; ++i) {
            std::cout << static_cast<char>('a' + m_rank[i][j]);
        }
        std::cout << "\n";
    }
}
