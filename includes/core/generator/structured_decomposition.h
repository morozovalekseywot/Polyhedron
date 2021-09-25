//
// Created by 159-mrv on 3/28/18.
//

#ifndef INFINIMESH_STRUCTURED_DECOMPOSITION_H
#define INFINIMESH_STRUCTURED_DECOMPOSITION_H

#include <allstd.h>

enum class Side;
class Configuration;

/**
 * @brief Функционал декомпозиции структурированной сетки
 *
 * Данный класс описывает "топологически" прямоугольную часть сетки (блок ячеек).
 * Для инициализации используется файл конфигурации и размеры глобальной сетки.
 *
 * Экземпляр класса содержит размеры блока, смещение блока в глобальной сетке,
 * число процессов по осям, а также координаты процесса.
 */
class StructuredDecomposition {
public:

    using UPtr = unique_ptr<StructuredDecomposition>;

    static StructuredDecomposition::UPtr create(const Configuration &config, array<size_t, 2> sizes) {
        return makeUnique<StructuredDecomposition>(config, sizes);
    };


    /** @brief Поддерживаемыые размерности декомпозиции */
    enum class Type { X, Y, XY };

    StructuredDecomposition() = default;

    /**
     * Определяет размерность декомпозиции по файлу конфигурации,
     * производит вычисление всех параметров декомпозиции.
     *
     * @param config Конфигурация
     * @param sizes Массив с количеством ячеек по размерностям [NX, NY, NZ] в глобальной сетке
     */
    explicit StructuredDecomposition(const Configuration &config, array<size_t, 2> sizes);


    /** @return Смещение блока по X (число столбцов ячеек левее данного процесса) */
    inline const size_t& x() const noexcept { return m_x; }

    /** @return Смещение блока по Y (число рядов ячеек ниже данного процесса) */
    inline const size_t& y() const noexcept { return m_y; }

    /** @return Размер блока по оси X */
    inline const size_t& nx() const noexcept { return m_nx; }

    /** @return Размер блока по оси Y */
    inline const size_t& ny() const noexcept { return m_ny; }


    /** @return Координата процесса по оси X */
    inline const int& pr_x() const noexcept { return m_pr_x; }

    /** @return Координата процесса по оси Y */
    inline const int& pr_y() const noexcept { return m_pr_y; }

    /** @return Число просессов по оси X */
    inline const int& pr_nx() const noexcept { return m_pr_nx; }

    /** @return Число процессов по оси Y */
    inline const int& pr_ny() const noexcept { return m_pr_ny; }


    /// @brief Ранг ячейки с указаными координатами
    /// @warning Ранг определяется корректно только для ячеек
    /// со смежных процессов!
    int cell_rank(size_t x, size_t y) const noexcept;


private:
    /**
     * @brief Находит в файле конфигурации тип декомпозици.
     * Не осуществляет проверок на возможность декомпозиции.
     * Инициализирует поле m_type.
     * @param config Указатель на конфигурацию
     */
    void read_config(const Configuration &config);


    /**
     * @brief Производит декомпозицию сетки по процессам,
     * определяет число процессов по каждой из осей,
     * а также координаты данного процесса.
     * Инициализирует поля m_pr_x, m_pr_y, m_pr_nx, m_pr_ny.
     * @param nx Число ячеек по оси X в глобальной сетке
     * @param ny Число ячеек по оси Y в глобальной сетке
     */
    void processes_decomposition(size_t nx, size_t ny);


    /** @brief Специализация для одномерной декомпозиции */
    void processes_decomposition_1d(int nx, int ny);


    /** @brief Специализация для двумерной декомпозиции */
    void processes_decomposition_2d(int nx, int ny);


    /**
     * @brief Заключительный этап декомпозиции.
     * Определяет размеры блока и смещение в глобальной сетке.
     * Инициализирует поля m_x, m_y, m_nx, m_ny.
     * @param nx Число ячеек по оси X в глобальной сетке
     * @param ny Число ячеек по оси Y в глобальной сетке
     */
    void mesh_decomposition(size_t nx, size_t ny);


    /** Тип разбиения */
    Type m_type;

    /** Координаты блока в глобальной сетке */
    size_t m_x, m_y;

    /** Размеры блока */
    size_t m_nx, m_ny;

    /** Координаты процесса */
    int m_pr_x, m_pr_y;

    /** Число процессов по осям */
    int m_pr_nx, m_pr_ny;
};

#endif //INFINIMESH_STRUCTURED_DECOMPOSITION_H
