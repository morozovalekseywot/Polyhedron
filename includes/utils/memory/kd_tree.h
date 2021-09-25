//
// Created by 159-egi on 2/2/18.
//

#ifndef INFINIMESH_KD_TREE_H
#define INFINIMESH_KD_TREE_H

#include <allstd.h>

#include "node_list.h"
#include "core/cell/cell.h"

class KDTree {
    using CList  = vector<Cell::Ptr>;
    using Ptr = shared_ptr<KDTree>;

    enum class SortBy {
        X, Y
    };

public:
    KDTree(CList nodes, int parts, int rank);
    KDTree(CList list, SortBy axis, int parts);

    void to_vector(vector<CList>& v);

private:

    SortBy opposite(SortBy axis);

    void sort(SortBy axis, CList& cells);

    CList m_list;

    vector<KDTree::Ptr> m_children;
};

#endif //INFINIMESH_KD_TREE_H
