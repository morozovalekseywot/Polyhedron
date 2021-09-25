//
// Created by 159-mrv on 1/17/19.
//

#include <utils/partitioning/voronoi_diagram.h>
#include <utils/math.h>

VoronoiDiagram::VoronoiDiagram(int size) {
    Vector2d zero = Vector2d::Zero(2);

    auto usize = (size_t)size;
    m_generators = vector<Vector2d>(usize, zero);
    m_adj_list = vector<vector<int>>(usize);
}

Vector2d VoronoiDiagram::get_generator(int index) const {
    return m_generators[index];
}

void VoronoiDiagram::set_generator(int index, Vector2d gen) {
    m_generators[index] = std::move(gen);
}

const vector<int>& VoronoiDiagram::adjacency(int index) const {
    return m_adj_list[index];
}

void VoronoiDiagram::update() {
    for(auto& adj: m_adj_list) {
        adj.clear();
    }

    int size = static_cast<int>(m_generators.size());


    for(int r = 0; r < size; ++r) {
        Vector2d v0 = m_generators[r];


        std::map<int, Vector2d> other_gens;
        for(int i = 0; i < size; ++i) {
            if (i != r) {
                other_gens[i] = m_generators[i];
            }
        }

        for(auto it = other_gens.begin(); it != other_gens.end(); ) {
            Vector2d kek = 0.5 * (v0 + it->second);
            double R = (kek - v0).squaredNorm();

            bool adj = true;
            for(auto it2 = other_gens.begin(); it2 != other_gens.end(); ++it2) {
                if (it == it2) {
                    continue;
                }

                double R2 = (it2->second - kek).squaredNorm();
                if (R2 < R) {
                    it = other_gens.erase(it);
                    adj = false;
                    break;
                }
            }

            if (adj) {
                ++it;
            }
        }

        for(auto& p: other_gens) {
            m_adj_list[r].push_back(p.first);
        }
    }
}

int VoronoiDiagram::cell_index(Vector2d vertex) const {
    double min = math::inf();
    int idx = -1;
    for(uint i = 0; i < m_generators.size(); ++i) {
        Vector2d diff = m_generators[i] - vertex;
        double dist = diff[0] * diff[0] + diff[1] * diff[1];

        if (dist < min) {
            min = dist;
            idx = i;
        }
    }
    return idx;
}
