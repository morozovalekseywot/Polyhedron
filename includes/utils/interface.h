//
// Created by 159-mrv on 7/17/19.
//

#ifndef INFINIMESH_INTERFACE_H
#define INFINIMESH_INTERFACE_H

#include <allstd.h>

/// @brief Структура для описания положения ампулы
class Interface {
public:
    using UPtr = unique_ptr<Interface>;

    inline static UPtr create(const vector<vector<Vector3d>>& points) {
        return makeUnique<Interface>(points);
    }

    Interface();

    ~Interface();

    explicit Interface(const vector<vector<Vector3d>>& points);


    void new_line_vertex(const Vector3d& v);

    void move(bool smooth = true);

    double sqr_distance(const Vector3d& v) const;


    void write(string filename) const;

    void reduce(const vector<Interface>& bords);


    inline const vector<vector<Vector3d>>& vertices() const { return m_points; }

    inline const Vector3d& vertex(size_t i, size_t j = 0) const { return m_points[i][j]; }

    inline double r_min() const { return m_r_min; }

    inline double r_max() const { return m_r_max; }

    inline double h2() const { return m_h2; }

private:

    // calc r_min, r_max, m_h2
    void calc_params();

    size_t m_size_x, m_size_y;
    vector<vector<Vector3d>> m_points;
    vector<vector<Vector3d>> m_offsets;
    vector<vector<double>> m_dists2;

    double m_r_min, m_r_max;
    double m_h2;
    double m_length;
};

#endif //INFINIMESH_INTERFACE_H
