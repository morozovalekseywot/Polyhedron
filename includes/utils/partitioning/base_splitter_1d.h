//
// Created by 159-mrv on 6/21/18.
//

#ifndef INFINIMESH_BASE_SPLITTER_1D_H
#define INFINIMESH_BASE_SPLITTER_1D_H

#include <allstd.h>
#include <utils/partitioning/splitter.h>

/// @brief Одномерный разделитель для декартовой/базовой сетки.
/// Для разбиения используется алгоритм бисекции.
class BaseSplitter1D : public Splitter {
public:

    /// @brief Конструктор класса
    /// @param size Размер одномерной сетки
    explicit BaseSplitter1D(size_t size);


    /// @brief Добавить нагрузку точке.
    /// @param point Координата точки
    /// @param weight Нагрузка
    void add_workload(size_t point, double weight);

    /// @brief Добавить нагрузку точке. Вариант, совместимый с базовым классом.
    /// @param point Целевая точка. Первая координата будет приведена к
    /// целочисленному типу, вторая будет проигнорирована.
    /// @param weight Нагрузка
    void add_workload(const SplitterFPoint& point, double weight) final;


    //// @brief Синхронизирует данные сплиттеров со всех процессов.
    void sync_data() final;


    /// @brief Производит разделение итнервала на несколько непересекающихся
    /// интервалов так, что нагрузки всех подинтервалов "примерно" совпадают.
    /// @param size Число подинтервалов
    void split(int size) final;


    /// @brief Вывести информацию о полученном разбиении: средняя нагрузка,
    /// дисбаланс и т. д.
    /// @param tabs
    void print_info(const string &tabs) const final;


    /// @brief Узнать какому процессу принадлежит точка.
    /// @param point Целевая точка.
    /// @return Ранг процесса.
    int rank(size_t point) const;


    /// @brief Узнать какому процессу принадлежит точка. Вариант, совместный
    /// с базовым классом.
    /// @param point Целевая точка. Будет использована только первая
    /// координата, приведенная к целочисленному типу.
    /// @return Ранг процесса.
    int rank(const SplitterFPoint& point) const final;


    /// @brief Вывод рангов в консоль (отладка).
    void print_ranks() const;


private:

    /// @brief Выполняет разбиение подобласти.
    /// @param b Прямоугольная подобласть
    /// @param r Диапозон рангов
    /// @param axis Приоритетная ось для разбиения (true -- X, false -- Y)
    void subsplit(size_t begin, size_t end, const RankRange &r);


    /// @brief Размеры сетки.
    size_t m_nx;

    /// @brief Количество подобластей для разбиения.
    unsigned int m_size;

    /// @brief Массив длины nx для хранения нагрузок.
    vector<double> m_workload;

    /// @brief Массив длины nx для хранения рангов.
    vector<int> m_rank;
};

#endif //INFINIMESH_BASE_SPLITTER_1D_H
