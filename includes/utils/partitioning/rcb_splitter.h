//
// Created by 159-mrv on 6/19/18.
//

#ifndef INFINIMESH_RCB_SPLITTER_H
#define INFINIMESH_RCB_SPLITTER_H

#include <allstd.h>
#include <utils/partitioning/splitter.h>

class BinaryTree  {
public:
    using UPtr = unique_ptr<BinaryTree>;

    BinaryTree(std::pair<SplitterFPoint, double>* points, size_t n_points, RankRange rank_range, uint axis);

    static BinaryTree::UPtr create(vector<std::pair<SplitterFPoint, double>>& points, int size) {
        return makeUnique<BinaryTree>(points.data(), points.size(), RankRange(0, size), 0);
    };

    int rank(const SplitterFPoint& point) const;

protected:
    uint m_coord_id;          // номер координаты, по которой делится область
    double m_coord;           // координата линии разбиения
    int m_rank;               // ранк процесса, если это лист разбиения

    BinaryTree::UPtr m_left;  // разбиение левого блока
    BinaryTree::UPtr m_right; // разбиение правого блока
};

/// @brief Реализация сплитера без ограничений.
class RCBSplitter : public Splitter {
public:

    RCBSplitter();

    /// @brief Добавить нагрузку точке. Если точка отсутствует, то она
    /// создается с указанной нагрузкой.
    /// @param point Координаты точки
    /// @param weight Нагрузка
    void add_workload(const SplitterFPoint& point, double weight) final;


    /// @brief Синхронизирует данные сплиттеров со всех процессов.
    void sync_data() final;


    /// @brief Производит разделение множества точек на несколько
    /// непересекающихся подмножеств так, что нагрузки всех подмножеств
    /// "примерно" совпадают.
    /// @param size Число множеств
    void split(int size) final;


    /// @brief Вывести информацию о полученном разбиении: средняя нагрузка,
    /// дисбаланс и т. д.
    /// @param tabs Отступ
    void print_info(const string &tabs) const final;


    /// @brief Узнать какому процессу принадлежит точка.
    /// @param point Целевая точка.
    /// @return Ранг процесса.
    int rank(const SplitterFPoint& point) const final;


protected:
    bool m_use_mpi;
    int m_size;
    unique_ptr<BinaryTree> m_tree;
    std::map<SplitterFPoint, double> m_points;
};

#endif //INFINIMESH_RCB_SPLITTER_H
