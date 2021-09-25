#pragma once

#include <utils/memory/node_list.h>

#include <problems/problem.h>

/// @brief Простой класс с искусственной задачей для проверки адаптации.
/// Задача заключается в адаптации сетки в соответсвии с заданной скалярной
/// функцией base_function (определена в классе). В результате адаптации
/// получаем уровни адаптации на сетке, равные целевой функции.
class AmrProblem : public Problem {

public:

    /// @brief Конструктор
    explicit AmrProblem(const Configuration& config);

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

    /// @brief Заполнить ячейки значениями целевой функции
    void stage_1(const NodeList::Part& cells);

    /// @brief Получить скалярное значение из ячейки
    double get_data(Cell_Ref cell) const;

    /// @brief Установить скалярное значение в ячейке
    void set_data(Cell_Ref cell, double value) const;

    /// @brief Целевая функция
    double base_function(const Vector3d& v, double t) const;

};
