//
// Created by 159-mrv on 8/27/19.
//

#ifndef INFINIMESH_CUBOID_GENERATOR_H
#define INFINIMESH_CUBOID_GENERATOR_H

#include <core/generator/structured_generator.h>

#if DIM3

class CuboidGenerator : public StructuredGenerator {

public:

    /// @brief Конструктор
    explicit CuboidGenerator(const Configuration &config);

    vector<vector<vector<Vertex_Ptr>>> create_vertices(Decomposition *decomp) const final;

    Cell_Ptr create_base_cell(size_t z) const final;

    Face_Ptr next_face_by_side(Cell_Ref cell, Side side) const final;

    void print_info(const string &tab) const final;


    double x_min() const;

    double y_min() const;

    double z_min() const;

    double x_max() const;

    double y_max() const;

    double z_max() const;

    double x_len() const;

    double y_len() const;

    double z_len() const;

private:
    double m_x_min, m_x_max, m_x_len;
    double m_y_min, m_y_max, m_y_len;
    double m_z_min, m_z_max, m_z_len;
};

#endif

#endif //INFINIMESH_CUBOID_GENERATOR_H
