//
// Created by 159-mrv on 8/14/18.
//

#ifndef INFINIMESH_GEOMETRY_H
#define INFINIMESH_GEOMETRY_H

#include <allstd.h>

class Cell;

/**************************************************************************//**
 * @brief Пространство имен для работы с геометрией, здесь собраны различные
 * функции для работы с различными геометрическими объектами и системами
 * координат.
 *****************************************************************************/
namespace geom {

    /// @brief Переводит пару декартовых координат в полярные
    /// @param vertex Декартовы координаты (x, y)
    /// @return Полярные координаты (r, phi)
    Vector2d to_polar(const Vector2d &vertex);

    /// @brief Переводит полярные координаты в декартовы
    Vector2d to_cartesian(double r, double phi);

    /// @brief Находит пересечение двух прямых
    /// @param v1, v2 - точки первой прямой
    /// @param v3, v4 - точки второй прямой
    /// @return Точка пересечения
    ///
    /// @warning Функция не осуществляет никаких проверок, если пересечение
    /// не единственно, то поведение не определено, также точка пересечения
    /// может не принадлежать отрезкам
    Vector2d lines_intersection(const Vector2d &v1, const Vector2d &v2,
                                const Vector2d &v3, const Vector2d &v4);

    /// @brief Квадрат расстояния от точки до отрезка
    /// Именно до отрезка! Если перпендикуляр, опущенный на прямую, находится
    /// вне отрезка, то выбирается концевая точка отрезка
    /// @param v Целевая точка
    /// @param v1, v2 Точки отрезка
    double sqr_distance(const Vector3d &v, const Vector3d &v1, const Vector3d &v2);

    /// @brief Квадрат расстояния от точки до площадки
    /// Именно до многоугольника! Если перпендикуляр, опущенный на плоскость,
    /// находится вне площадки, то выбирается граничная точка
    /// @param v Целевая точка
    /// @param v1, v2, v3, v4 Вершины четырехугольника
    double sqr_distance(const Vector3d &v, const Vector3d &v1, const Vector3d &v2,
                        const Vector3d &v3, const Vector3d &v4);

    /// @brief Вычисляет площадь полигона
    /// @param vertices Точки, упорядоченные против часовой стрелки
    /// @return Площадь полигона
    /// @warning Не рекомендуется в приложениях с параллельностью по тредам
    /// из-за использования динамического массива vector
    double polygon_area(const vector<Vector2d>& vertices);

    /// @brief Вычисляет площадь треугольника
    /// @param vertices Вершины треугольника
    double polygon_area(const array<Vector2d, 3>& vertices);

    /// @brief Вычисляет площадь треугольника
    /// @param vertices Вершины треугольника
    double polygon_area(const array<Vector3d, 3>& vertices);

    /// @brief Вычисляет площадь четырехугольника
    /// @param vertices Точки, упорядоченные против часовой стрелки
    double polygon_area(const array<Vector2d, 4>& vertices);

    /// @brief Вычисляет площадь четырехугольника
    /// @param vertices Точки, упорядоченные против часовой стрелки
    double polygon_area(const array<Vector3d, 4>& vertices);

    /// @brief Внешняя нормаль четырехугольного полигона
    Vector3d normal(const array<Vector3d, 4>& vertices);

    /// @brief Находит точку, симметричную относительно прямой
    /// @param p Заданая точка
    /// @param v1, v2 Точки, задающие прямую
    Vector2d symmetric_point(const Vector2d& p, const Vector2d& v1, const Vector2d& v2);

    /// @brief Разделяет секторальную ячейку на две части прямой {v1, v2}.
    /// Возвращает объемную долю части ячейки, которая находится выше прямой.
    /// "выше/ниже" расчитывается по оси радиуса: точка плоскости находится выше
    /// линии {v1, v2}, если отрезок, соединяющий её с началом координат
    /// пересекает линию {v1, v2}.
    /// Если вся ячейка выше прямой, тогда возвращается 1.0, если ниже - 0.0.
    /// Таким образом, возвращаемое значения всегда лежит на отрезке [0.0, 1.0]
    ///
    /// @param cell Целевая секторальная ячейка
    /// @param v1, v2 Точки прямой
    /// @return Объемная часть ячейки, находящаяся выше разделяющей линии
    double split_cell(shared_ptr<Cell> cell, const Vector2d& v1, const Vector2d& v2);
}

#endif //INFINIMESH_GEOMETRY_H
