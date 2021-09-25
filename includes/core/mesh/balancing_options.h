//
// Created by 159-mrv on 10/15/18.
//

#ifndef INFINIMESH_BALANCING_OPTIONS_H
#define INFINIMESH_BALANCING_OPTIONS_H

#include <allstd.h>

class Configuration;

/// @brief Типы балансировки.
enum class BalancingType : int {
    NONE, /// @brief Балансировка отсутствует.
    X,    /// @brief Разбиение вдоль оси X.
    Y,    /// @brief Разбиение вдоль оси Y.
    XY    /// @brief Двумерное разбиение.
};

/// @brief Наборы ячеек, по которым проводится разбиение.
enum class CellsSet : int {
    BASE, /// @brief Только базовые ячейки.
    LEAF  /// @brief Листовые ячейки.
};

/// @brief Опции балансировки нагрузки.
struct BalancingOptions {

    /// @brief Конструктор по умолчанию.
    BalancingOptions() = default;

    /// @brief Конструктор класса, читает конфигурацию и инициализирует все
    /// параметры.
    explicit BalancingOptions(const Configuration& config);

    inline bool use_balancing() const noexcept {
        return type != BalancingType::NONE;
    }

    /// @brief Тип балансировки.
    BalancingType type;

    /// @brief Ячейки, по которым проводится балансировка.
    CellsSet cells;

    /// @brief Число процессов по оси X при ОДНОМЕРНОМ разбиении.
    int proc_per_x;

    /// @brief Число процессов по оси Y при ОДНОМЕРНОМ разбиении.
    int proc_per_y;
};

#endif //INFINIMESH_BALANCING_OPTIONS_H
