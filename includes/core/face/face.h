#pragma once

#include <allstd.h>
#include <core/face/face_flag.h>
#include <core/cell/side.h>
#include <core/indexed_object.h>

class Vertex;
class Cell;

/**************************************************************************//**
 * @brief Класс для хранения ориентированной грани.
 *
 * @details Поддерживается как хранение одномерной, так и двумерной
 * граней. Информация о размерности грани заключена в количестве точек.
 * Так, одномерная грань содержит две точки, двумерная - четыре.
 *
 * @details Грань является ориентированной, а следовательно имеет две нормали.
 * outward_normal - Внешняя нормаль, в двумерном случае является правой
 * нормалью к вектору грани, в трехмерном случае определяется по правилу
 * правой руки при обходе четырех вершин грани.
 * inward_normal  - Внутренняя нормаль, противоположна внешней.
 *
 * @details Помимо граней хранит слабые указатели на соседние ячейки.
 * outward/inward_cell - соседи грани в направлении соответствующих нормалей.
 *****************************************************************************/
class Face : public std::enable_shared_from_this<Face>,
             public IndexedObject
             {

    /// @{
    /// @brief Сокращения, используемые в описании класса.
    using Vertex_Ptr = shared_ptr<Vertex>;
    using Vertex_Ref = const shared_ptr<Vertex> &;
    using Cell_Ptr   = shared_ptr<Cell>;
    using Cell_Ref   = const shared_ptr<Cell> &;
    using Cell_WPtr   = weak_ptr<Cell>;
    /// @}

public:
    /// @brief Сокращенное именование shared_ptr на грань.
    using Ptr = shared_ptr<Face>;

    /// @brief Константаная ссылка на shared_ptr для передачи в функции.
    using Ref = const shared_ptr<Face> &;



    /** @{ ****************************************************************//**
     *
     * @name Конструторы грани
     *
     * @brief Конструкторы принимают по два параметра
     * @param begin Указатель на начало
     * @param end Указатель на конец

     * @warning Порядок вершин имеет значение, поскольку грани являются
     * ориентированными, объяснение в описании, см. Face.
     *
     * @throw runtime_error, если одна из вершин nullptr
     *
     *************************************************************************/

    /// @brief Конструктор одномерной грани
    Face(Vertex_Ref v1, Vertex_Ref v2);

    /// @brief Создает новую грань и возвращает умный указатель на неё
    static inline Face::Ptr create(Vertex_Ref v1, Vertex_Ref v2) {
        return make_shared<Face>(v1, v2);
    }

    /// @brief Конструктор двумерной грани
    Face(Vertex_Ref v1, Vertex_Ref v2, Vertex_Ref v3, Vertex_Ref v4);

    /// @brief Создает новую грань и возвращает умный указатель на неё
    static inline Face::Ptr create(Vertex_Ref v1, Vertex_Ref v2, Vertex_Ref v3, Vertex_Ref v4) {
        return make_shared<Face>(v1, v2, v3, v4);
    }

    /// @brief Переворачивает уже созданную сторону: меняет порядок точек,
    /// направление нормалей, inward и outward ячейки
    void reverse();

#if DIM3
    /// @brief Поворачивает грань вокруг внешней нормали, актуально только для
    /// двумерных граней
    /// @param r Количество смещений -4 < r < 4
    void rotate(int r);
#endif

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Характеристики грани.
     *
     *************************************************************************/

    /// @brief Площадь грани (длина в 2D)
    double area() const;

    /// @brief Площадь для осесимметричных задач
    double area_as() const;

    /// @brief Центр грани (= среднее вершин)
    const Vector3d &center() const;

    /// @brief Центр грани
    Vector2d center_2d() const;

    /// @brief Нормаль к грани, проведенная в сторону соседа
    /// @param neighbor Соседняя ячейка
    /// @throw runtime_error, если грань не имеет такого соседа
    const Vector3d &normal(Cell_Ref neighbor) const;

    /// @brief Нормаль к грани, проведенная в сторону соседа
    /// @param neighbor Соседняя ячейка
    /// @throw runtime_error, если грань не имеет такого соседа
    const Vector3d &normal(Cell* neighbor) const;

    /// @brief Нормаль к грани, проведенная в сторону соседа
    /// @param neighbor Соседняя ячейка
    /// @throw runtime_error, если грань не имеет такого соседа
    Vector2d normal_2d(Cell_Ref neighbor) const;

    /// @brief Внешняя нормаль к грани
    const Vector3d &outward_normal() const;

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Работа с флагом граничных условий
     *
     *************************************************************************/

    /// @brief Получить флаг граничных условий
    FaceFlag flag() const;

    /// @brief Устанавить флаг граничных условий
    void set_flag(FaceFlag flag);

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Get/Set методы для вершин грани.
     *
     * По мере возможностей рекомендуется использовать простые указатели.
     *
     *************************************************************************/

    /// @brief Возвращает ссылку на массив вершин
    const vector<Vertex_Ptr>& vertices() const noexcept;

    /// @brief Число вершин грани
    size_t size() const noexcept;

    /// @brief Получить вершину
    Vertex_Ref vertex(int idx) const noexcept;

    /// @brief Установить новую вершину одномерной грани
    void set_vertex(int idx, Vertex_Ref v) noexcept;

    /// @brief Сырой указатель на вершину
    Vertex* vertex_raw(int idx) const noexcept;

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Get/Set методы для прилегающих ячеек.
     *
     *************************************************************************/

    /// @brief Внешняя ячейка (по нормали)
    Cell_Ptr inward_cell() const;

    /// @brief Установить внешнюю ячейку
    void set_inward_cell(Cell_Ref cell);

    /// @brief Внутренняя ячейка
    Cell_Ptr outward_cell() const;

    /// @brief Уставновить внутренюю ячейку
    void set_outward_cell(Cell_Ref cell);

    /// @brief Существует ли у cell сосед через данную грань.
    bool neighbor_exist(const Cell *cell) const noexcept;

    /// @brief Существует ли у cell сосед через данную грань.
    bool neighbor_exist(Cell_Ref cell) const noexcept;

    /// @brief Возвращает соседа cell через данную грань.
    /// @throw runtime_error, если ячейка cell не смежна грани или сосед не
    /// существует.
    Cell* neighbor_raw(const Cell *cell) const;

    /// @brief Возвращает соседа cell через данную грань.
    /// @throw runtime_error, если ячейка cell не смежна грани или сосед не
    /// существует.
    Cell_Ptr neighbor(const Cell *cell) const;

    /// @brief Возвращает соседа cell через данную грань.
    /// @throw runtime_error, если ячейка cell не смежна грани или сосед не
    /// существует.
    Cell_Ptr neighbor(Cell_Ref cell) const;

    /** @} */


private:

    /// @brief Инициализация одномерной грани. Вызывается в конструктуре,
    /// вычисляет площадь, нормаль и центр грани.
    void initialize();

    vector<Vertex_Ptr> m_vertices;

    FaceFlag m_flag;

    Cell_WPtr m_inward_cell;
    Cell_WPtr m_outward_cell;

    double   m_area;
    Vector3d m_center;
    Vector3d m_inward_normal;
    Vector3d m_outward_normal;

};

/// @brief Строковое представление грани (выводит вершины)
string to_string(const Face * face);

/// @brief Строковое представление грани (выводит вершины)
string to_string(Face::Ref face);


/// @brief Структура грани для сериализации.
struct ClearFace {
    IndexedObject::index_type index;
    IndexedObject::index_type vertex[4];
    IndexedObject::index_type inward_cell;
    IndexedObject::index_type outward_cell;
    int flag;

    ClearFace() = default;

    explicit ClearFace(Face::Ref face);

    explicit ClearFace(const string &line);
};

/// @brief Вывод грани в ASCII поток
ostream& operator<<(ostream& out, const ClearFace& face);
