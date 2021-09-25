#pragma once

#include <allstd.h>

class Face;
class Face;
class Cell;
class Mesh;
class Configuration;

enum class AdaptationFlag : int;

/// @brief Интерфейс задачи
class Problem {
public:
    using Face_Ptr = shared_ptr<Face>;
    using Face_Ref = const shared_ptr<Face>&;
    using Cell_Ptr = std::shared_ptr<Cell>;
    using Cell_Ref = const std::shared_ptr<Cell>&;

public:
    using UPtr = unique_ptr<Problem>;

    Problem();

    static Problem::UPtr create(const Configuration &config);


    /// @return Размер данных ячейки, определяется задачей
    virtual uint cell_data_size() const = 0;


    /// @brief Ининциализация ячеек сетки начальными значениями
    /// @param mesh Сетка
    virtual void initialization(Mesh *mesh) = 0;


    /// @brief Главная функция решателя
    /// @param mesh Сетка
    /// @return Чистое время решения, без учета обменов.
    virtual double solution_step(Mesh *mesh) = 0;


    /// @brief Возвращает внутрирасчетное время
    double current_time() const;


    /// @brief Установить внутрирасчетное время
    void set_time(double time);



    /** @{ ****************************************************************//**
     *
     * @name Методы, связаные с адаптацией.
     *
     *************************************************************************/

    /// @brief Определяет, необходимо ли ячейке адаптироваться
    /// @param cell Целевая ячейка
    /// @return Флаг адаптации: COARSE, NONE или SPLIT
    virtual AdaptationFlag adaptation_criterion(const shared_ptr<Cell>& cell) = 0;

    /// @brief Распределяет данные родительской ячейки между детьми при разбиении
    virtual void split_data(const shared_ptr<Cell>& parent, const vector<shared_ptr<Cell>> &children) = 0;

    /// @brief Вычисляет данные родительской ячейки при огрублении
    virtual void coarse_data(const shared_ptr<Cell>& parent, const vector<shared_ptr<Cell>> &children) = 0;

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Методы вывода информации.
     *
     *************************************************************************/

    /// @brief Печатает некоторые данные в консоль
    /// @param tab Отступ перед выходными данными
    /// @param real_time Внутрирассчетное время
    virtual void print_info(const char *tab) const = 0;

    /// @brief Позволяет получить из ячейки параметр для записи
    /// @param cell Целевая ячейка
    /// @return Значение параметра
    virtual double get_cell_param(const shared_ptr<Cell>& cell, const string& name) const = 0;

    /// @brief Позволяет получить параметры солвера для записи
    /// @return Параметры солвера на запись, порядок соответствует выдачи params_to_write()
    virtual double get_integral_param(const string& name) const = 0;

    /** @} */

protected:

    /// @brief Внутрирасчетное время
    double m_time;
};
