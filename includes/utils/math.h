//
// Created by egi on 12/16/17.
//

#ifndef INFINIMESH_MATH_H
#define INFINIMESH_MATH_H

#include <allstd.h>

/// @brief В данном пространстве имен собраны некоторые необходимые
/// математические функции.
namespace math {

    /// @brief Целочисленный логарифм
    size_t int_log(double x, double base);

    /// @brief Компоратор для std::pair<size_t, size_t>
    struct pairhash {
    public:
        // Магия по ссылке:
        // https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
        size_t operator()(const std::pair<size_t, size_t> &x) const {
            auto lhs = std::hash<size_t>{}(x.first);
            auto rhs = std::hash<size_t>{}(x.second);

            lhs ^= rhs + 0x9e3779b97f4a7c16 + (lhs << 6) + (lhs >> 2);

            return lhs;
        }
    };


    /// @brief Возвращает положительную бесконечность.
    inline double inf() noexcept {
        return std::numeric_limits<double>::infinity();
    }


    /// @brief Возвращает среднее арифметическое.
    double mean(const vector<double>& values);


    /// @brief Возвращает среднеквадратическое отклонение.
    /// @param values Значения.
    /// @param mean Среднее арифметическое (если посчитано).
    double std_deviation(const vector<double>& values, double mean = NAN);


    /// @brief Степень двойки, 2^n
    template <class T>
    inline T pow2(T n) { return 1 << n; }


    /// @brief Степень четверки, 4^n
    template <class T>
    inline T pow4(T n) { return 1 << (n << 1); }
}

#endif // INFINIMESH_MATH_H
