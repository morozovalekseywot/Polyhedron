//
// Created by 159-mrv on 5/7/18.
//

#ifndef INFINIMESH_BASE_SPLITTER_2D_H
#define INFINIMESH_BASE_SPLITTER_2D_H

#include <allstd.h>
#include "splitter.h"

/// @brief Двумерный разделитель для декартовой/базовой сетки.
/// Для разбиения используется алгоритм recursive coordinate bisection.
class BaseSplitter2D : public Splitter {
public:

    /// @brief Конструктор класса
    /// @param sizes Размеры сетки
    explicit BaseSplitter2D(const SplitterIPoint& sizes);


    /// @brief Добавить нагрузку точке.
    /// @param point Координаты точки
    /// @param weight Нагрузка
    void add_workload(const SplitterIPoint& point, double weight);


    /// @brief Добавить нагрузку точке. Вариант, совместимый с базовым классом.
    /// @param point Целевая точка. Первая координата будет приведена к
    /// целочисленному типу, вторая будет проигнорирована.
    /// @param weight Нагрузка
    void add_workload(const SplitterFPoint& point, double weight) final;


    /// @brief Заполнение случайными нагрузками (отладка)
    void random_workload();


    /// @brief Синхронизирует данные сплиттеров со всех процессов.
    void sync_data();


    /// @brief Производит разделение множества точек на несколько
    /// непересекающихся подмножеств так, что нагрузки всех подмножеств
    /// "примерно" совпадают.
    /// @param size Число множеств/
    void split(int size);


    /// Вывести информацию о полученном разбиении: средняя нагрузка,
    /// дисбаланс и т. д.
    /// @param tabs
    void print_info(const string &tabs) const;


    /// Узнать какому процессу принадлежит точка.
    /// @param point Целевая точка.
    /// @return Ранг процесса.
    int rank(const SplitterIPoint& point) const;


    /// @brief Узнать какому процессу принадлежит точка. Вариант, совместный
    /// с базовым классом.
    /// @param point Целевая точка. Будет использована только первая
    /// координата, приведенная к целочисленному типу.
    /// @return Ранг процесса.
    int rank(const SplitterFPoint& point) const final;


    /// @brief Вывод рангов в консоль (отладка).
    void print_ranks() const;


private:

    /// @brief Структура для хранения двумерного прямоугольного блока.
    struct Block {
        size_t x1, y1;
        size_t x2, y2;

        Block(size_t _x1, size_t _y1, size_t _x2, size_t _y2)
                : x1(std::min(_x1, _x2)),
                  y1(std::min(_y1, _y2)),
                  x2(std::max(_x1, _x2)),
                  y2(std::max(_y1, _y2)) {

        }

        size_t x_size() const {
            return x2 - x1;
        }

        size_t y_size() const {
            return y2 - y1;
        }
    };


    /// @brief Выполняет разбиение подобласти.
    /// @param b Прямоугольная подобласть
    /// @param r Диапозон рангов
    /// @param axis Приоритетная ось для разбиения (true -- X, false -- Y)
    void subsplit(const Block& b, const RankRange& r, bool axis = true);


    /// @brief Размеры сетки.
    size_t m_nx, m_ny;

    /// @brief Количество подобластей для разбиения.
    unsigned int m_size;

    /// @brief Массив (nx, ny) для хранения нагрузок.
    vector<vector<double>> m_workload;

    /// @brief Массив (nx, ny) для хранения рангов.
    vector<vector<int>> m_rank;
};

#endif //INFINIMESH_BASE_SPLITTER_2D_H
