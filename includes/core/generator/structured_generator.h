//
// Created by 159-mrv on 3/28/18.
//

#ifndef INFINIMESH_STRUCTURED_GENERATOR_H
#define INFINIMESH_STRUCTURED_GENERATOR_H

#include <allstd.h>
#include <core/cell/side.h>
#include <core/generator/structured_decomposition.h>
#include <core/generator/geometry_generator.h>


/**************************************************************************//**
 *
 * @brief Сеточный генератор для двумерной структурированной сетки
 *
 *****************************************************************************/
class StructuredGenerator : public GeometryGenerator {

public:

    /** @{ ****************************************************************//**
     *
     * @brief Функции свойственные непосредственно генератору
     * структурированной сетки
     *
     *************************************************************************/

    /// @brief Конструктор по умолчанию, ничего не заполняет
    /// (а все потому что m_global_nx или m_global_ny могут вычисляться
    /// алгоритмически в зависимости от реализации)
    explicit StructuredGenerator(GeometryType type);


    /// @brief Возвращает массив базовых ячеек сетки.
    /// Возвращаются только ячейки, принадлежащие данному процессу согласно
    /// проведенной декомпозиции.
    /// Возвращаемые ячейки связаны друг с другом и проиндексированы по
    /// кривой, заполняющей пространство.
    /// Для корректной работы необходима реализация функции create_vertices
    vector<Cell_Ptr> create_cells(Decomposition *decomp) const override;

    /// @brief Возвращает двумерный массив вершин
    virtual vector<vector<vector<Vertex_Ptr>>> create_vertices(Decomposition *decomp) const = 0;


    /// @brief Число ячеек вдоль оси X
    size_t global_nx() const noexcept;

    /// @brief Число ячеек вдоль оси Y
    size_t global_ny() const noexcept;

    /// @brief Число ячеек вдоль оси Z
    size_t global_nz() const noexcept;

    /// @brief Полное число ячеек
    size_t n_cells() const noexcept final;


    /// @brief Возвращает положение соседа ячейки со стороны side
    /// на кривой, заполняющей пространство
    size_t neighbor_z(Cell_Ref cell, Side side) const final;


    /// @brief Переводит координаты декартовой сетки в уникальную координату
    size_t xyz2u(size_t x, size_t y, size_t z = 0) const;

    /// @brief Переводит уникальную координату в декартовы
    array<size_t, 2> u2xy(size_t u) const;

    /// @brief Переводит уникальную координату в декартовы
    array<size_t, 3> u2xyz(size_t u) const;

    /** @} */

protected:
    size_t m_global_nx;
    size_t m_global_ny;
    size_t m_global_nz;
};

#endif //INFINIMESH_STRUCTURED_GENERATOR_H
