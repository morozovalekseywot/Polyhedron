#pragma once

#include <utils/memory/node_list.h>

#include <problems/problem.h>

struct CellData {
public:
    CellData() = default;

    double vol; ///< Объемная доля
    std::array<double, 3> n; ///< Внешняя нормаль
};


class SurfaceProblem : public Problem {
public:

    /// @brief Конструктор
    explicit SurfaceProblem(const Configuration& config);

    /// @brief Размер данных для хранения в ячейке
    uint cell_data_size() const final;

    /// @brief Заполнение сетки начальными данными
    void initialization(Mesh *mesh) final;

    /// @brief Шаг решателя
    double solution_step(Mesh *mesh) final;



    /// @brief Критерий адаптации, возращает флаг адаптации в зависимости от
    /// состояния в ячейке
    AdaptationFlag adaptation_criterion(Cell_Ref cell) final;

    /// @brief Оператор переноса данных при разбиении ячейки
    void split_data(Cell_Ref parent, const vector<Cell_Ptr> &children) final;

    /// @brief Оператор переноса данных при огрублении ячейки
    void coarse_data(Cell_Ref parent, const vector<Cell_Ptr> &children) final;



    /// @brief Вывести информацию в консоль после шага
    void print_info(const char *tab) const final;

    /// @brief Значение параметра с именем name в ячейке
    double get_cell_param(Cell_Ref cell, const string& name) const final;

    /// @brief Значение интегрального параметра с именем name
    double get_integral_param(const std::string& name) const final;

private:

    CellData get_state(Cell_Ref cell) const;
    void set_state(Cell_Ref cell, const CellData& value) const;
};
