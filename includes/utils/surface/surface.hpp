/// @file Я опишу структуру всех классов в одном файле, но классы необходимо
/// будет разделить по разным файлам, а также разделить объявление и
/// определение функций в заголовочные .h файлы и .cpp файлы.
/// Любая функция располагается внутри некоторого класса.
/// Классы именуются с использованием CamelCase, функции именуются
/// с использованием snake_case. Определения из enum class
/// прописными буквами. Члены класса начинаются с префикса m_ за
/// исключением полей небольших структур к которым предполагается
/// публичный доступ.

// не забываем include guarding
#pragma once

#include <vector>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <utils/surface/triangle.hpp>
#include <utils/surface/face.hpp>
#include <utils/surface/bbox.hpp>
#include <boost/container/static_vector.hpp>
#include <Dense> /// from Eigen library
#include <algorithm>


using boost::container::static_vector;
using Eigen::Vector3d;
using Eigen::Matrix3d;

/// using для остальных классов снаружи запрещен!


/// @brief Используем namespace для всех классов, поскольку типы вроде
/// Vertex и Triangle могут встречаться в других местах проекта
namespace surf
{
    /// @brief Класс поверхности, версия для овощей
    /// Есть только конструктор и функция is_inside, реализованная
    /// первоклассниками. В целом, этого достаточно для решения задачи.
    class Surface
    {
    #define equal_d(a, b) abs(a - b) < m_length / 1e8

    public:

        Surface() : m_length(0.0) {}
        /// @brief Конструктор класса из stl файла
        /// @param filename Имя stl файла
        explicit Surface(const std::string &filename);

        /// @return True, если точка v находится внутри поверхности
        bool is_inside(const Vector3d &v) const;

        /// @brief Центрировать (центр BoundingBox в начале координат)
        void centering();

        /// @brief Поворот вокруг оси на угол
        void rotate(Vector3d n, double phi);

        /// Геттеры сделаны для тестов, можно потом удалить
        double getMLength() const;

        const std::vector<Vertex> &getMVertices() const;

        const std::vector<Triangle> &getMTriangles() const;

    private:
        double m_length;                    ///< Характерный размер
        std::vector<Vertex> m_vertices;     ///< Список вершин
        std::vector<Triangle> m_triangles;  ///< Список треугольников
    };

}