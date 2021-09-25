//
// Created by 159-mrv on 6/9/18.
//

#ifndef INFINIMESH_SPLITTER_H
#define INFINIMESH_SPLITTER_H

#include <allstd.h>

/// @brief Простая структура для хранения диапозона рангов.
/// Используется классический полуинтервал [begin, end).
struct RankRange {
    int begin, end;

    RankRange(int r1, int r2) noexcept :
            begin(r1), end(r2) { }

    inline int size() const noexcept { return end - begin; }
};


/// @brief Типы точек.
using SplitterIPoint = array<size_t, 2>;
using SplitterFPoint = array<double, 2>;


/// @brief Функция выводит статистку по загрузке процессов.
void print_workload_info(const string &tabs, const vector<double>& workloads);


/// @brief Полностью абстрактный splitter.
/// Разбивает множество точек с нагрузками на компактные непересекающиеся
/// подмножества таким образом, что нагрузки всех подмножеств "примерно"
/// совпадают.
class Splitter {

public:

    /// @brief Сокращенное именование указателя.
    using UPtr = unique_ptr<Splitter>;


    /// @brief Добавить нагрузку точке.
    /// Если точка отсутствует, то она создается с указанной нагрузкой.
    /// @param point Координаты точки
    /// @param weight Нагрузка
    virtual void add_workload(const SplitterFPoint& point, double weight) = 0;


    /// @brief Синхронизирует данные сплиттеров со всех процессов.
    virtual void sync_data() = 0;


    /// @brief Производит разделение множества точек на несколько
    /// непересекающихся подмножеств так, что нагрузки всех подмножеств
    /// "примерно" совпадают.
    /// @param size Число множеств
    virtual void split(int size) = 0;


    /// @brief Вывести информацию о полученном разбиении: средняя нагрузка,
    /// дисбаланс и т. д.
    /// @param tabs Отступ
    virtual void print_info(const string &tabs) const = 0;


    /// @brief Узнать какому процессу принадлежит точка.
    /// @param point Целевая точка.
    /// @return Ранг процесса.
    virtual int rank(const SplitterFPoint& point) const = 0;
};

#endif //INFINIMESH_SPLITTER_H
