//
// Created by 159-mrv on 9/6/19.
//

#ifndef INFINIMESH_SECTOR3D_GENERATOR_H
#define INFINIMESH_SECTOR3D_GENERATOR_H

#include <core/generator/structured_generator.h>

#if DIM3

class Sector3dGenerator : public StructuredGenerator {
public:
    explicit Sector3dGenerator(const Configuration &config);

    vector<vector<vector<Vertex_Ptr>>> create_vertices(Decomposition *decomp) const final;

    Cell_Ptr create_base_cell(size_t z) const final;

    Face_Ptr next_face_by_side(Cell_Ref cell, Side side) const final;

    void print_info(const string &tab) const final;

    double r_in() const;

    double r_out() const;

    double angle1() const;

    double angle2() const;

private:
    double m_global_r1;
    double m_global_r2;
    double m_alpha;
    double m_beta;
};

#endif

#endif //INFINIMESH_SECTOR3D_GENERATOR_H
