//
// Created by egi on 12/24/17.
//

#ifndef INFINIMESH_LAYERCONTROL_H
#define INFINIMESH_LAYERCONTROL_H

#include <allstd.h>

enum TAGS {
    COUNT_TAG, DATA_TAG, ADAPTATION_FLAGS_TAG, CELLS_INDICES_TAG, PRIMITIVES_COUNT_TAG, PRIMITIVES_INDICES_TAG, DOUBLE_TAG
};

class Cell;
class NodeList;

class LayerControl {
public:
    explicit LayerControl(shared_ptr<NodeList> layer);

    virtual ~LayerControl() = default;

    shared_ptr<NodeList> cells();

    void rebuild();

    int rank();

protected:

    shared_ptr<NodeList> m_layer;
    vector<shared_ptr<Cell>> m_sorted_cells;
};

#endif //INFINIMESH_LAYERCONTROL_H
