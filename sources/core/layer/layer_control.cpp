//
// Created by egi on 12/24/17.
//

#include <core/cell/cell.h>
#include <core/layer/layer_control.h>
#include <utils/memory/node_list.h>

LayerControl::LayerControl(NodeList::Ptr layer)
    : m_layer(std::move(layer)) {

}

NodeList::Ptr LayerControl::cells() {
    return m_layer;
}

int LayerControl::rank() {
    return m_layer->rank();
}

void LayerControl::rebuild() {
    m_sorted_cells.clear();

    for (auto cell: *m_layer) {
        m_sorted_cells.push_back(cell);
    }

    std::sort(m_sorted_cells.begin(), m_sorted_cells.end(),
              [](Cell::Ref cell1, Cell::Ref cell2) -> bool {
                  return cell1->id() < cell2->id();
              });
}