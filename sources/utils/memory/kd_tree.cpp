//
// Created by 159-egi on 2/2/18.
//

// C++
#include <algorithm>
#include <cmath>
#include <iostream>

// InfiniMesh
#include "utils/memory/kd_tree.h"
#include "core/vertex/vertex.h"
#include "core/cell/data_holder.h"


KDTree::KDTree(CList cells, int parts, int rank)
{
    if(parts > 1) {
        auto init_sort = SortBy::X;

        sort(init_sort, cells);

        auto part_left  = (int)std::round(parts / 2.0);
        int part_right = parts - part_left;

        auto half_size = (size_t)std::round(part_left * cells.size() / (double)parts);

        m_children.push_back(make_shared<KDTree>(CList(cells.begin(), cells.begin() + half_size), opposite(init_sort), part_left));
        m_children.push_back(make_shared<KDTree>(CList(cells.begin() + half_size, cells.end()), opposite(init_sort), part_right));
    }
    else {
        m_list = std::move(cells);
    }
}

KDTree::KDTree(CList list, SortBy axis, int parts)
{
    if(parts > 1) {
        sort(axis, list);

        auto part_left  = (int)std::round(parts / 2.0);
        int part_right = parts - part_left;

        auto half_size = (size_t)std::round(part_left * list.size() / (double)parts);

        m_children.push_back(make_shared<KDTree>(CList(list.begin(), list.begin() + half_size), opposite(axis), part_left));
        m_children.push_back(make_shared<KDTree>(CList(list.begin() + half_size, list.end()), opposite(axis), part_right));
    }
    else {
        m_list = std::move(list);
    }
}

KDTree::SortBy KDTree::opposite(SortBy axis) {
    if(axis == SortBy::X) { return SortBy::Y; }
    else                  { return SortBy::X; }
}

void KDTree::sort(SortBy axis, CList& cells) {
    if(axis == SortBy::X) {
        std::sort(cells.begin(), cells.end(), [](const Cell::Ptr& a, const Cell::Ptr& b) -> bool {
            return a->center_2d()[0] > b->center_2d()[0];
        });
    }
    else if(axis == SortBy::Y){
        std::sort(cells.begin(), cells.end(), [](const Cell::Ptr& a, const Cell::Ptr& b) -> bool {
            return a->center_2d()[1] > b->center_2d()[1];
        });
    }
    else {
        throw runtime_error("Error! Can't sort by this dimension");
    }
}

void KDTree::to_vector(vector<CList>& v) {
    if(m_children.empty()) {
        v.push_back(m_list);
    }
    else {
        for(auto& child: m_children) {
            child->to_vector(v);
        }
    }
}

