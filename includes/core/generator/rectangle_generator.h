//
// Created by 159-mrv on 8/16/18.
//

#ifndef INFINIMESH_RECTANGLE_GENERATOR_H
#define INFINIMESH_RECTANGLE_GENERATOR_H

#include <core/generator/structured_generator.h>

#define ROTATION_ANGLE (0.0*M_PI/4.0)

class RectangleGenerator : public StructuredGenerator {

public:

    /// @brief Конструктор
    explicit RectangleGenerator(const Configuration &config);

    vector<vector<vector<Vertex_Ptr>>> create_vertices(Decomposition *decomp) const final;

    Cell_Ptr create_base_cell(size_t z) const final;

    Face_Ptr next_face_by_side(Cell_Ref cell, Side side) const final;

    void print_info(const string &tab) const final;


    double x_min() const;

    double y_min() const;

    double x_max() const;

    double y_max() const;

    double x_len() const;

    double y_len() const;

private:
    double m_x_min, m_x_max;
    double m_y_min, m_y_max;
    double m_x_len, m_y_len;
};

#endif //INFINIMESH_RECTANGLE_GENERATOR_H
