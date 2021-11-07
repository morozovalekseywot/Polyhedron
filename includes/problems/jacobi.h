#pragma once

#include <utils/memory/node_list.h>
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

// Форма нашего объекта
struct Figure
{
    std::vector<std::pair<Vector3d, Vector3d>> sides; // пара: (нормаль, центр)

    Vector3d right_up = {0.7, 0.2, 0.0};
    Vector3d right_down = {0.7, -0.2, 0.0};
    Vector3d left_up = {-0.70, 0.2, 0.0};
    Vector3d left_down = {-0.7, 0.2, 0.0};

    Figure(const Vector3d &r_up = {0.7, 0.2, 0.0}, const Vector3d &r_down = {0.7, -0.2, 0.0},
           const Vector3d &l_down = {-0.7, 0.2, 0.0}, const Vector3d &l_up = {-0.70, 0.2, 0.0}) : right_up(r_up), right_down(r_down),
                                                                                                  left_down(l_down), left_up(l_up)
    {
        Vector3d up = right_up - left_up;
        Vector3d right = right_down - right_up;
        Vector3d down = left_down - right_down;
        Vector3d left = left_up - left_down;

        sides.resize(4);
        sides[0] = std::make_pair(Vector3d{up.x(), up.y(), 0.0}, up / 2);
        sides[1] = std::make_pair(Vector3d{right.x(), right.y(), 0.0}, right / 2);
        sides[2] = std::make_pair(Vector3d{down.x(), down.y(), 0.0}, down / 2);
        sides[3] = std::make_pair(Vector3d{left.x(), left.y(), 0.0}, left / 2);
        for (auto &side: sides)
        {
            if (side.first.x() * side.second.x() + side.first.y() * side.second.y() +
                side.first.z() * side.second.z() < 0.0)
            {
                side.first = -side.first;
            }
        }
    }

    double is_inside(const Vector3d &r)
    {
        for (auto &side: sides)
        {
            Vector3d vec = side.second - r;
            if (side.first.x() * vec.x() + side.first.y() * vec.y() +
                side.first.z() * vec.z() < 0.0)
                return 0.0;
        }

        return 1.0;
    }

    ~Figure() = default;
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

    JacobiCellData get_state(Cell_Ref cell) const;

    void set_state(Cell_Ref cell, const JacobiCellData &value) const;

    Figure figure;
};