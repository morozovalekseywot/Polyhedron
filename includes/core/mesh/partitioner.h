//
// Created by 159-mrv on 11/1/18.
//

#ifndef INFINIMESH_PARTITIONER_H
#define INFINIMESH_PARTITIONER_H

#include <allstd.h>

#include <utils/partitioning/splitter.h>
#include <core/mesh/balancing_options.h>


class Cell;
class Splitter;
class StructuredGenerator;


class Partitioner {
private:
    using Cell_Ref = const shared_ptr<Cell>&;

public:
    Partitioner(const BalancingOptions& options, StructuredGenerator* geometry);

    void add_workload(Cell_Ref cell);

    void sync_data();

    void split();

    void print_info(const string& tab) const;

    int rank(Cell_Ref cell) const;

private:

    BalancingOptions m_options;

    StructuredGenerator* m_geometry;

    unique_ptr<Splitter> m_splitter;

    std::function<SplitterFPoint(Cell_Ref)> m_location;

    std::function<double(Cell_Ref)> m_workload;
};

#endif //INFINIMESH_PARTITIONER_H
