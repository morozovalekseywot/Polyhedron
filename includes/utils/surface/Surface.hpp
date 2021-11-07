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
#include <utils/surface/Triangle.hpp>
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
// TODO разделить по файлам
    class Surface
    {
    public:

        /// @brief Конструктор класса из stl файла
        /// @param filename Имя stl файла
        Surface(const std::string &filename)
        {
            /// В конструкторе необходимо заполнить массивы m_vertices
            /// и m_triangles.

            /// Считываем треугольники из файла в промежуточный массив,
            /// где каждый треугольник представляется набором из трех
            /// вершин v1, v2, v3 и нормали n, далее работаем с ним,
            /// файл закрываем.
            string path = "../examples";
            path = path + filename;
            ifstream file(path, std::ios_base::in);
            if (!file.is_open()) // если файл не был открыт
            {
                throw runtime_error("file isn't open\n");
            }
            std::vector<Face> faces;
//            std::vector<std::array<Vector3d, 4>> raw_triangles;
            while (!file.eof())
            {
                string t;
                getline(file, t);
                getline(file, t);
                Vector3d a, b, c;
                file >> a.x() >> a.y() >> a.z();
                file >> b.x() >> b.y() >> b.z();
                file >> c.x() >> c.y() >> c.z();
                faces.emplace_back(Face({a, b, c}));
                getline(file, t);
                getline(file, t);
                getline(file, t);
            }
            file.close();

            /// По ходу считывания из файла вычисляем m_length. Более точно:
            /// строим BoundingBox, затем m_length считаем как длину диагонали.
            /// В дальнейшем m_length это характерный размер.
            BBox bbox(faces);
            m_length = bbox.len.norm();
            /// Далее самая сложная операция, удалить дубликаты вершин.

            /// Вариант 1. Тупой, за O(N^2).
            /// Проходим по треугольникам в raw_triangles, треугольник кидаем
            /// в m_triangles, затем у каждого треугольника проходим по вершинам.
            /// Сначала ищем вершину в m_vertices (простым перебором), если
            /// вершины там нет (с точностью до eps*m_length), тогда добавляем
            /// новую вершину в конец m_vertices, а индекс этой вершины
            /// записываем в вершину к m_triangles.back().

            /// Вариант 2. Адекватный, за O(N log(N)).
            /// Вершины добавляются параллельно в m_vertices и в некоторый
            /// set, для которого пишется определенный компаратор.
            /// Вместо прохода по всем добавленным вершинам, новая вершина
            /// просто ищется в set.
            /// Компаратор следующий: вершины эквивалентны, если совпадают
            /// с точностью до eps*m_length, иначе сравниваются их проекции
            /// на некоторую довольно случайную ось.
            /// Хотя у меня есть сомнения, что вы осознаете, что тут написано.
            /// Пусть eps = 1e-8

            auto comp = [this](const Eigen::Vector3d &a, const Eigen::Vector3d &b) -> bool
            {
                if (equal_d(a.x(), b.x()))
                {
                    if (equal_d(a.y(), b.y()))
                    {
                        if (equal_d(a.z(), b.z()))
                            return false;
                        else
                            return a.z() < b.z();
                    } else
                        return a.y() < b.y();
                } else
                    return a.x() < b.x();
            };

            std::set<Eigen::Vector3d, decltype(comp)> set_vtx(comp);
            for (auto &face: faces)
            {
                set_vtx.insert(face.vects[0]);
                set_vtx.insert(face.vects[1]);
                set_vtx.insert(face.vects[2]);
            }

            m_vertices.resize(set_vtx.size());
            int size = 0;
            for (auto &v: set_vtx)
                m_vertices[size++] = Vertex(v);

            if (!std::is_sorted(m_vertices.begin(), m_vertices.end(), comp))
            {
                throw std::runtime_error("Write to Lesha"); // TODO delete this
            }

            size = 0;
            m_triangles.resize(faces.size());
            for (auto &face: faces)
            {
                // ищем индексы наших вершин в массиве m_vertices, сверху всталвена проверка на то, что
                // отсортрован, если там упадёт, то нельзая бин поиском искать (lower_bound)
                auto it1 = std::lower_bound(m_vertices.begin(), m_vertices.end(), face.vects[0], comp);
                auto it2 = std::lower_bound(m_vertices.begin(), m_vertices.end(), face.vects[1], comp);
                auto it3 = std::lower_bound(m_vertices.begin(), m_vertices.end(), face.vects[2], comp);
                size_t idx1 = std::distance(m_vertices.begin(), it1);
                size_t idx2 = std::distance(m_vertices.begin(), it2);
                size_t idx3 = std::distance(m_vertices.begin(), it3);
                m_triangles[size] = Triangle(size, {idx1, idx2, idx3}, face.normal);
                size++;
            }
            /// Ладно, пофиг, реализовали первый варинат. Теперь у нас есть
            /// m_vertices только с положениями, без дубликатов и
            /// m_triangles, которые на них ссылаются.

            // TODO посмотрите на метод set_next_triangle и доделайте конструктор

            /// Проходим по треугольникам в m_triangles, у них по вершинам,
            /// для вершин записываем смежные треугольники.

            /// Проходим по треугольникам в m_triangles, проходим по их
            /// вершинам, у вершин в треугольниках ищем смежные треугольники.
        }

        /// @return True, если точка v находится внутри поверхности
        bool is_inside(const Vector3d &v) const
        {
            /// Здесь простая реализация. Будем использовать простейшее правило
            /// четный-нечетный (even-odd rule). Из какой-нибудь точки выпускается
            /// луч в произвольном направлении. Затем отыскиваются пересечения луча
            /// и треугольников поверхности. Если пересечений нечетное количество,
            /// значит точка лежит внутри поверхности, иначе -- снаружи.
            /// Если точка и луч достаточно случайные, то можно считать, что
            /// вероятность попадания луча на границу треугольника равна нулю
            /// (все должно работать короче).
        }

    private:
        double m_length;                    ///< Характерный размер
        std::vector<Vertex> m_vertices;     ///< Список вершин
        std::vector<Triangle> m_triangles;  ///< Список треугольников
    };

}