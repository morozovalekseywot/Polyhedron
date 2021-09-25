//
// Created by 159-mrv on 12/6/18.
//

#ifndef INFINIMESH_REAL_SECTOR_GENERATOR_H
#define INFINIMESH_REAL_SECTOR_GENERATOR_H

#include <core/generator/geometry_generator.h>


class RealSectorGenerator : public GeometryGenerator {
    using Face_Ptr = shared_ptr<Face>;

    struct Index {
        int part;
        size_t x;
        size_t y;
    };

public:
    explicit RealSectorGenerator(const Configuration &config);

    vector<Cell_Ptr> create_cells(Decomposition *decomp) const final;

    Cell_Ptr create_base_cell(size_t z) const final;

    Face_Ptr next_face_by_side(Cell_Ref cell, Side side) const final;

    /// @brief Полное число ячеек
    size_t n_cells() const noexcept final;

    /// @brief Возвращает положение соседа ячейки со стороны side
    /// на кривой, заполняющей пространство
    size_t neighbor_z(Cell_Ref cell, Side side) const final;

    /// @brief Возвращает положение соседа ячейки со стороны side
    Index neighbor(const Index& index, Side side) const;

    void print_info(const string &tab) const final;


    double r_min() const noexcept {
        return m_r1;
    }

    double r_out() const noexcept {
        return m_r2;
    }

    double angle() const noexcept {
        return m_alpha;
    }

private:

    size_t center_z(size_t x_base, size_t y_base) const;

    size_t sector_z(size_t x_base, size_t y_base) const;

    size_t get_z(const Index& index) const;


    void create_vertices();

    Face_Ptr get_center_hface(size_t x, size_t y) const;

    Face_Ptr get_center_vface(size_t x, size_t y) const;

    Face_Ptr get_sector_hface(size_t x, size_t y) const;

    Face_Ptr get_sector_vface(size_t x, size_t y) const;




    /// @brief Раствор угла
    double m_alpha;

    /// @brief Радиус начала структурированной части
    double m_r1;

    /// @brief Внешний радиус
    double m_r2;



    const double m_limiter = 3 * M_PI_4;

    /// @brief Число больших неструктурированных квадратов, если угол раствора
    /// меньше m_limiter, тогда равно 1, иначе 2
    size_t m_nc;

    /// @brief Разбиение больших неструктурированных квадратов ( = Nx / (2 nc) )
    size_t m_nx;

    /// @brief Разбиение по углу в структурированной части
    size_t m_Nx;

    /// @brief Разбиение по радиусу в структурированной части
    size_t m_Ny;

    /// @brief Полное число ячеек ( = nc * nx * (nx + 2 * Ny ) )
    size_t m_size;


    vector<vector<Vertex_Ptr>> m_center_vertices;
    vector<vector<Vertex_Ptr>> m_sector_vertices;
};


#endif //INFINIMESH_REAL_SECTOR_GENERATOR_H
