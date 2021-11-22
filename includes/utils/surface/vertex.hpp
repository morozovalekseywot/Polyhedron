#pragma once

#include <vector>
#include <string>
#include <map>
#include <set>

#include <boost/container/static_vector.hpp>
#include <Dense> /// from Eigen library

using boost::container::static_vector;
using Eigen::Vector3d;

// TODO убрать все using ... из Vertex.hpp, Surface.hpp и остальных файлов, которые я добавил (ctrl+R - поиск с заменой)
namespace surf
{
    /// @class Описание вершины. Хранит положение и список смежных треугольников
    class Vertex
    {
    public:

        Vertex() = default;

        explicit Vertex(const Vector3d &point);

        ~Vertex() = default;

        /// @brief Максимальное количество треугольников, разделяющих данную
        /// вершину, если вы не найдете адекватные файлы с моделями, то есть
        /// возникнет необходимость использовать большое значение max_triangles
        /// (> 20), тогда придется использовать vector вместо static_vector
        static const int max_triangles = 10;

        /// Метод чтобы в конце конструктора расставить для вершин смежные треуголники
        void set_next_triangle(size_t idx);

        /// конвертирование в Eigen::Vector3d
        operator Eigen::Vector3d() const;

        // Возможно эти операторы вообще не нужны
        Vector3d operator+(const Vertex &a) const;

        Vector3d operator-(const Vertex &a) const;

        double mod() const;

        void normalize();

        template<class T>
        Vector3d operator*(const T &a) const;

        template<class T>
        Vector3d operator/(const T &a) const;

        double operator*(const Vertex &a) const;

        Vector3d operator-();


        /// @brief Массив с индексами треугольников, которые разделяют данную
        /// вершину
        /// @details Контейнер static_vector имитирует поведение std::vector за
        /// исключением того, что может быть увеличен только до max_triangles.
        /// Благодаря известному размеру, поле типа static_vector располагается
        /// в памяти последовательно с другими полями, а не в куче (куча это
        /// термин, см. heap) как std::vector. Короче по различным причинам
        /// это будет работать быстрее для небольших значений max_triangles.
        static_vector<size_t, max_triangles> triangles;

        /// @brief Непосредственное положение вершины
        Vector3d v;
    };
}

/// Векторное произведение
inline Vector3d cross(const surf::Vertex &a, const surf::Vertex &b);

inline std::ostream &operator<<(std::ostream &os, const surf::Vertex &a);