//
// Created by 159-mrv on 1/17/19.
//

#ifndef INFINIMESH_VORONOI_DIAGRAM_H
#define INFINIMESH_VORONOI_DIAGRAM_H

#include <allstd.h>

class VoronoiDiagram {

public:

    explicit VoronoiDiagram(int size);

    Vector2d get_generator(int index) const;

    void set_generator(int index, Vector2d gen);

    const vector<int>& adjacency(int index) const;

    /// Бичарский алгоритм, работает за O(N^3), где N - число точек.
    void update();


    int cell_index(Vector2d vertex) const;

private:

    vector<Vector2d> m_generators;

    vector<vector<int>> m_adj_list;

};


#endif //INFINIMESH_VORONOI_DIAGRAM_H
