//
// Created by 159-egi on 2/20/17.
//

#ifndef INFINIMESH_CALCULATION_CONTROLLER_H
#define INFINIMESH_CALCULATION_CONTROLLER_H

#include <allstd.h>

class MeshWriter;
class Mesh;
class Problem;

/**
 * @brief Менеджер процесса расчёта
 * @details Управляет счётом, определяет периодичность сохранения, необходимость адаптации.
 * Содержит основной цикл работы приложения.
 */
class CalculationController
{
public:
    CalculationController(const string& config_filename, const string& checkpoint_filename);
    ~CalculationController();

    /// @brief Запускает основной цикл работы приложения
    void run();

private:

    double solve();

    void perf_test();

    /// @brief Печатает информацию о текущем шаге расчёта
    void print_step_info(double step_time, long full_time);

    /// @brief Выводит статистику о времени рассчета
    void print_time_stat(long calc_time, long write_time, long solve_time,
                         long adapt_time, long balance_time);

    /// @brief Флаги команд времени выполнения.
    enum class CommandFlag {
        EXIT, ///< @brief Преждевременное завершение.
        NONE  ///< @brief Нет побочных действий.
    };

    /// @brief Выполняет команду времени выполнения.
    CommandFlag exec_command();


    /// @brief Внутрирасчётное время, при достижении которого выполняется остановка.
    double         m_max_time;

    /// @brief Шаг расчёта при достижении которого выполняется остановка.
    size_t         m_max_step;

    /// @brief Шаг расчёта.
    size_t         m_step;

    /// @brief Флаг, необходима ли адаптация?
    bool           m_adaptation;

    /// @brief Частота адаптации.
    size_t         m_adapt_frequency;

    /// @brief Частота вызова балансировки.
    uint           m_balancing_frequency;

    string m_command_filename;          /// @brief Файл с командами.
    vector<double> m_checkpoints;       /// @brief Время записи чекпоинта.
    string         m_checkpoint_prefix; /// @brief Название файла сохранения.
    bool           m_checkpoint_binary; /// @brief Использовать бинарный формат.

    unique_ptr<Problem>    m_problem;   /// @brief Решатель.
    unique_ptr<Mesh>       m_mesh;      /// @brief Сетка.
    std::vector<std::unique_ptr<MeshWriter>> m_writers;    /// @brief Сохраняет сетку.
};

#endif //INFINIMESH_CALCULATION_CONTROLLER_H
