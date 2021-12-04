#pragma once

#include <utils/memory/node_list.h>
#include <utils/surface/surface.hpp>
#include <core/face/face.h>
#include <problems/problem.h>

// В каждой ячейке будем хранить индекс, значения скалярного потенциала
// на двух шагах u1 и u2 (чтобы выводить итерации методов), компоненты скорости Vx,Vy, Vz.
struct JacobiCellData
{
    double idx; ///< индекс
    double u1; ///< значение скалярного потенциала на предыдущем шаге
    double u2; ///< значение скалярного потенциала на следующем шаге
    Vector3d v; ///< компоненты скорости
    double vol; ///< Объёмная доля

    JacobiCellData() : idx(0), u1(0.0), u2(0.0), v({0.0, 0.0, 0.0}), vol(0.0)
    {}

    JacobiCellData(double idx, double u1, double u2, const Vector3d &v, double vol) : idx(idx), u1(u1), u2(u2), v(v), vol(vol)
    {}

    JacobiCellData(const JacobiCellData &data) : idx(data.idx), u1(data.u1),
                                                 u2(data.u2), v(data.v), vol(data.vol)
    {}

    ~JacobiCellData() = default;
};

class Jacobi : public Problem
{
public:
    /// @brief Конструктор
    explicit Jacobi(const Configuration &config);

    /// @brief Размер данных для хранения в ячейке
    uint cell_data_size() const override;

    /// @brief Заполнение сетки начальными данными
    void initialization(Mesh *mesh) final;

    /// @brief Шаг решателя
    double solution_step(Mesh *mesh) override;

    /// @brief Критерий адаптации, возращает флаг адаптации в зависимости от
    /// состояния в ячейке
    AdaptationFlag adaptation_criterion(const shared_ptr<Cell> &cell) override;

    /// @brief Оператор переноса данных при разбиении ячейки
    void split_data(const shared_ptr<Cell> &parent, const vector<shared_ptr<Cell>> &children) override;

    /// @brief Оператор переноса данных при огрублении ячейки
    void coarse_data(const shared_ptr<Cell> &parent, const vector<shared_ptr<Cell>> &children) override;

    /// @brief Вывести информацию в консоль после шага
    void print_info(const char *tab) const override;

    /// @brief Значение параметра с именем name в ячейке
    double get_cell_param(const shared_ptr<Cell> &cell, const string &name) const override;

    /// @brief Значение интегрального параметра с именем name
    double get_integral_param(const string &name) const override;

private:
    double m_eps = 0.0; // max|u2-u1|
    double m_delta = 0.0;
    bool first_step = true; // был ли совершён первый шаг расчёта

    JacobiCellData get_state(Cell_Ref cell) const;

    void set_state(Cell_Ref cell, const JacobiCellData &value) const;

    surf::Surface figure;
};