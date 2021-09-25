//
// Created by 159-mrv on 6/21/18.
//

#include <control/mpi_wrapper.h>

#include <utils/partitioning/base_splitter_1d.h>


BaseSplitter1D::BaseSplitter1D(size_t size)
        : m_nx(size), m_size(0) {
    m_workload = vector<double>(m_nx, 0.0);
    m_rank = vector<int>(m_nx, -1);
}

void BaseSplitter1D::add_workload(size_t point, double weight) {
    m_workload[point] += weight;
}

void BaseSplitter1D::add_workload(const SplitterFPoint& point, double weight) {
    m_workload[static_cast<size_t>(point[0])] += weight;
}

void BaseSplitter1D::sync_data() {
    m_workload = mpi::sum(m_workload);
}

void BaseSplitter1D::split(int size) {
    m_size = static_cast<uint>(size);
    RankRange r(0, m_size);
    subsplit(0, m_nx, r);
}

void BaseSplitter1D::subsplit(size_t begin, size_t end, const RankRange &r) {
    int n_rank = r.size();

    //mpi::cout << "subsplit: " << begin << " " << end << " " << r.begin << " " << r.end << "\n";

    if (n_rank == 1) {
        for (size_t i = begin; i < end; ++i)
            m_rank[i] = r.begin;
        return;
    }

    vector<double> wl_sum(end - begin, 0.0);

    wl_sum[0] = m_workload[begin];
    for (size_t i = 1; i < wl_sum.size(); ++i) {
        wl_sum[i] = wl_sum[i - 1] + m_workload[begin + i];
    }

    int min_range = n_rank / 2;
    int max_range = n_rank - min_range;

    double min_diff = 1.0 / 0.0;
    size_t delemeter = wl_sum.size() * min_range / max_range;
    int left_rank_range = min_range;

    for (size_t i = 1; i < wl_sum.size(); ++i) {
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

        //mpi::cout << "  " << i << " " << m_workload[begin + i] << " "
        //          << wl_sum[i] << " " << diff1 << " " << diff2 << " " << min_diff << "\n";
    }

    //mpi::cout << "  del: " << delemeter << " " << left_rank_range << "\n";


    RankRange left_range(r.begin, r.begin + left_rank_range);
    subsplit(begin, begin + delemeter, left_range);

    RankRange right_range(r.begin + left_rank_range, r.end);
    subsplit(begin + delemeter, end, right_range);
}

int BaseSplitter1D::rank(size_t point) const {
    return m_rank[point];
}

int BaseSplitter1D::rank(const SplitterFPoint &point) const {
    return m_rank[static_cast<size_t>(point[0])];
}

void BaseSplitter1D::print_ranks() const {
    for (size_t j = 0; j < m_nx; ++j) {
        std::cout << static_cast<char>('a' + m_rank[j]);
    }
    std::cout << std::endl;
}

void BaseSplitter1D::print_info(const string &tabs) const {
    vector<double> wl(m_size, 0.0);
    for (size_t i = 0; i < m_nx; ++i) {
        wl[m_rank[i]] += m_workload[i];
    }
    print_workload_info(tabs, wl);
}