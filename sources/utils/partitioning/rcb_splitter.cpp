//
// Created by 159-mrv on 6/19/18.
//

#include <control/mpi_wrapper.h>
#include <utils/partitioning/rcb_splitter.h>

RCBSplitter::RCBSplitter() :
        m_use_mpi(false),
        m_size(mpi::size()) {

}

void RCBSplitter::add_workload(const SplitterFPoint& point, double weight) {
    m_points[point] += weight;
}

void RCBSplitter::sync_data() {
    m_use_mpi = true;
}

// Так, я знаю, что это дно вариант, но щито поделать,
// возможно, когда-нибудь дойдут руки сделать нормально.
void RCBSplitter::split(int size) {
    m_size = size;

    vector<std::pair<SplitterFPoint, double>> points;
    if (!m_use_mpi) {
        // Splitter используется на одном процессе
        // или на всех процессах независимо
        points.reserve(m_points.size());
        for (const auto &p: m_points) {
            points.emplace_back(p);
        }
    } else {
        // Здесь я собираю точки со всех процессов в данный
        vector<double> data(3 * m_points.size());
        vector<double> full_data = mpi::all_gather(data);
        data.clear();

        points.reserve(full_data.size() / 3);
        for (size_t i = 0; i < full_data.size(); i += 3) {
            points.emplace_back(
                    std::make_pair(
                            SplitterFPoint({full_data[i], full_data[i + 1]}),
                            full_data[i + 2]
                    ));
        }
        full_data.clear();
    }

    m_tree = BinaryTree::create(points, m_size);
}

int RCBSplitter::rank(const SplitterFPoint& point) const {
    return m_tree->rank(point);
}


int BinaryTree::rank(const SplitterFPoint& point) const {
    if (m_left && m_right) {
        return (point[m_coord_id] < m_coord) ?
               m_left->rank(point) : m_right->rank(point);
    } else {
        return m_rank;
    }
};

BinaryTree::BinaryTree(std::pair<SplitterFPoint, double> *points, size_t n_points, RankRange rank_range, uint axis)
        : m_coord_id(axis), m_coord(0.0), m_rank(0),
          m_left(nullptr), m_right(nullptr) {

    // База рекурсии, нечего делить,
    // всё уже проинициализировано
    if (rank_range.size() == 1) {
        m_coord_id = axis;
        m_rank = rank_range.begin;
        return;
    }

    using namespace std::chrono;

    auto begin = steady_clock::now();

    int n_ranks = rank_range.size();

    std::sort(points, points + n_points, [axis](const std::pair<SplitterFPoint, double> &a,
                                                const std::pair<SplitterFPoint, double> &b) -> bool {
        return a.first[axis] < b.first[axis];
    });


    vector<double> wl_sum(n_points, 0.0);
    wl_sum[0] = points[0].second;

    for (size_t i = 1; i < n_points; ++i) {
        wl_sum[i] = wl_sum[i - 1] + points[i].second;
    }

    int min_range = n_ranks / 2;
    int max_range = n_ranks - min_range;


    double min_diff = 1.0 / 0.0;
    size_t delemeter = n_points * min_range / max_range;
    int left_rank_range = min_range;

    auto sigma = 1.0e-4 * fabs(points[n_points].first[axis] - points[0].first[axis]);
    for (size_t i = 1; i < n_points; ++i) {
        if (fabs(points[i - 1].first[axis] - points[i].first[axis]) < sigma)
            continue;

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

    uint next_axes = (axis + 1) % 2;

    m_coord = points[delemeter].first[axis];
    m_rank = 0;
    m_left = makeUnique<BinaryTree>(points, delemeter,
                                     RankRange(rank_range.begin, rank_range.begin + left_rank_range),
                                     next_axes);


    m_right = makeUnique<BinaryTree>(points + delemeter, n_points - delemeter,
                                      RankRange(rank_range.begin + left_rank_range, rank_range.end),
                                      next_axes);

    auto end = steady_clock::now();
    std::cout << "n_points: " << n_points << ". time: " << 1000.0 * duration<double>(end - begin).count() << " ms\n";
}


void RCBSplitter::print_info(const string &tabs) const {
    vector<double> wl(m_size, 0.0);
    for(auto &p: m_points) {
        wl[rank(p.first)] += p.second;
    }

    print_workload_info(tabs, wl);
}