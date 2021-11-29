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

/// using для остальных классов снаружи запрещен!


/// @brief Используем namespace для всех классов, поскольку типы вроде
/// Vertex и Triangle могут встречаться в других местах проекта
namespace surf
{
#define equal_d(a, b) abs(a - b) < m_length / 1e8

    /// @brief Класс поверхности, версия для овощей
    /// Есть только конструктор и функция is_inside, реализованная
    /// первоклассниками. В целом, этого достаточно для решения задачи.
    class Surface
    {
    public:

        /// @brief Конструктор класса из stl файла
        /// @param filename Имя stl файла
        explicit Surface(const std::string &filename);

        /// @return True, если точка v находится внутри поверхности
        bool is_inside(const Vector3d &v) const
        {
            /// Здесь простая реализация. Будем использовать простейшее правило
            /// четный-нечетный (even-odd rule). Из какой-нибудь точки выпускается
            /// луч в произвольном направлении. Затем отыскиваются пересечения луча
            /// и треугольников поверхности.
//            Vector3d s = m_vertices[rand() % m_vertices.size()].v; // конец луча
            Vector3d s = {0.0, 0.0, 0.0};
            for (auto &idx: m_triangles[m_vertices.back().triangles[0]].vertices) // берём треугольник смежный последней точке
                s = s + m_vertices[idx].v;
            s /= 3.0;
            Vector3d l = s - v; // направляющий вектор прямой
            s = v + 1.1 * l; // продлил луч
            l = s-v;

            int count = 0; // количество пересечений

            auto scalar = [](const Vector3d &a, const Vector3d &b) -> double
            {
                return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
            };

            for (auto &triangle: m_triangles)
            {
//                В первую очередь, проверяется, пересекает ли контрольный луч плоскость просматриваемого
//                треугольника: точки v и 𝑆 должны лежать с разных сторон или скалярные произведения
//                нормали с радиус-векторами от любой вершины до этих точек должны иметь разные знаки. В
//                противном случае – переход к следующему треугольнику без дальнейших проверок.
                Eigen::Vector3d first = m_vertices[triangle.vertices[0]].v; // TODO починить
                if (scalar(triangle.normal, v - first) * scalar(triangle.normal, s - first) > 0.0)
                    continue;

//                ур-е плоскости: Ax + By + Cz + D = 0
//                ур-е прямой: x = v.x + m * t, y = v.y + p * t, z = v.y + l * t
//                A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2)
//                B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2)
//                C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2)
//                -D = x1 (y2 z3 - y3 z2) + x2 (y3 z1 - y1 z3) + x3 (y1 z2 - y2 z1)
                std::array<Vector3d, 3> arr;
                for (int i = 0; i < 3; i++)
                    arr[i] = m_vertices[triangle.vertices[i]].v;

                double A, B, C, D;
                A = arr[0].y() * (arr[1].z() - arr[2].z()) + arr[1].y() * (arr[2].z() - arr[0].z()) + arr[2].y() * (arr[0].z() - arr[1].z());
                B = arr[0].z() * (arr[1].x() - arr[2].x()) + arr[1].z() * (arr[2].x() - arr[0].x()) + arr[2].z() * (arr[0].x() - arr[1].x());
                C = arr[0].x() * (arr[1].y() - arr[2].y()) + arr[1].x() * (arr[2].y() - arr[0].y()) + arr[2].x() * (arr[0].y() - arr[1].y());
//                D = -(arr[0].x() * (arr[1].y() * arr[2].z() - arr[1].z() * arr[2].y()) + arr[1].x() * (arr[2].y() * arr[0].z() - arr[2].z() * arr[0].y()) +
//                      arr[2].x() * (arr[0].y() * arr[1].z() - arr[0].z() * arr[1].y()));
                D = -A * arr[0].x() - B * arr[0].y() - C * arr[0].z();

                for (auto &p: arr)
                    if (abs(A * p.x() + B * p.y() + C * p.z() + D) > 1e-5)
                        throw std::runtime_error("wrong plane equation"); // TODO убрать проверку, если всё работает

                double t = -(D + A * v.x() + B * v.y() + C * v.z()) / (A * l.x() + B * l.y() + C * l.z());
                if (t < 0.0)
                    continue;
                Vector3d x = {v.x() + t * l.x(), v.y() + t * l.y(), v.z() + t * l.z()};
                double area_ = ((x - arr[0]).cross(x - arr[1])).norm() + ((x - arr[0]).cross(x - arr[2])).norm() + ((x - arr[1]).cross(x - arr[2])).norm();
                if (abs(area_ - 2 * triangle.area) > triangle.area * 1e-5)
                    continue; // точка не в треугольнике
                else
                    count++;
            }

            /// Если пересечений нечетное количество,
            /// значит точка лежит внутри поверхности, иначе -- снаружи.
            /// Если точка и луч достаточно случайные, то можно считать, что
            /// вероятность попадания луча на границу треугольника равна нулю
            /// (все должно работать короче).

            return count % 2 != 0;
        }

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