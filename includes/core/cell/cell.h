#pragma once

#include <allstd.h>

#include <core/indexed_object.h>
#include <utils/memory/node.h>

//#define USE_POLAR_MESH

class Vertex;
class Face;
class Face;
class FacesList;
class DataHolder;
enum class Side;
enum class GeometryType : int;

enum class AdaptationFlag : int { COARSE = -1, NONE = 0, SPLIT = 1 };

/**************************************************************************//**
 * @brief Класс ячейки.
 *
 * Необходимо максимально упростить.
 * Весь сложный функционал должен быть вынесен в управляющую структуру,
 * то есть класс Mesh.
 * Класс Cell должен предоставлять только простейшие get/set методы
 * и некоторое число проверок
 *****************************************************************************/
class Cell : public std::enable_shared_from_this<Cell>,
             public IndexedObject,
             public Node
{

    /// @{
    /// @brief Сокращения, используемые в определении класса
    using Vertex_Ptr = shared_ptr<Vertex>;
    using Vertex_Ref = const shared_ptr<Vertex> &;
    using Face_Ptr   = shared_ptr<Face>;
    using Face_Ref   = const shared_ptr<Face> &;
    using Faces_Ptr  = shared_ptr<FacesList>;
    using Faces_Ref  = const shared_ptr<FacesList> &;
    using Data_Ptr   = shared_ptr<DataHolder>;
    using Data_Ref   = const shared_ptr<DataHolder> &;
    /// @}

public:
    using IDType = size_t;

    /** @{ ****************************************************************//**
     * @name Сокращенные именования указателей.
     *************************************************************************/
    using Ptr   = shared_ptr<Cell>;
    using Ref   = const shared_ptr<Cell> &;

    using WPtr  = weak_ptr<Cell>;
    using WRef  = const weak_ptr<Cell> &;

    using UPtr  = unique_ptr<Cell>;
    using URef  = const unique_ptr<Cell> &;
    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Конструкторы.
     *
     * @brief Все конструкторы принимают в качестве аргументов указатели на
     * грани (или списки граней FacesList).
     *
     * Аргументы инициализируют поля m_left_faces, m_right_faces, и т.д.,
     * рассчитываются характеристики ячейки (центр, объем и др),
     * а все остальные поля класса выставляются значениями по умолчанию.
     * Связи между ячейкой и гранями также остаются не выставленными.
     *
     *************************************************************************/

    Cell(Face_Ref left, Face_Ref bottom, Face_Ref right, Face_Ref top);

    Cell(Faces_Ref left, Faces_Ref bottom, Faces_Ref right, Faces_Ref top);

    static inline Ptr create(Face_Ref left, Face_Ref bottom, Face_Ref right, Face_Ref top) {
        return make_shared<Cell>(left, bottom, right, top);
    }

    static inline Ptr create(Faces_Ref left, Faces_Ref bottom, Faces_Ref right, Faces_Ref top) {
        return make_shared<Cell>(left, bottom, right, top);
    }

    Cell(Face_Ref left, Face_Ref bottom, Face_Ref right,
         Face_Ref top, Face_Ref front, Face_Ref back);

    Cell(Faces_Ref left, Faces_Ref bottom, Faces_Ref right,
         Faces_Ref top, Faces_Ref front, Faces_Ref back);

    static inline Ptr create(Face_Ref left, Face_Ref bottom, Face_Ref right,
                             Face_Ref top, Face_Ref front, Face_Ref back) {
        return make_shared<Cell>(left, bottom, right, top, front, back);
    }

    static inline Ptr create(Faces_Ref left, Faces_Ref bottom, Faces_Ref right,
                             Faces_Ref top, Faces_Ref front, Faces_Ref back) {
        return make_shared<Cell>(left, bottom, right, top, front, back);
    }

    ~Cell() final = default;

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Примитивы
     *
     * @brief Get/Set методы для работы с примитивами ячейки: гранями,
     * вершинами, а также соседями.
     *
     *************************************************************************/

    /// @brief Правильно ли направлена грань на стороне side, грань считается
    /// правильно ориентированной, если левая нормаль грани является внешней
    /// к данной ячейке
    bool is_correct_face(Side side) const;

    bool is_correct_face(Face_Ref face) const;

    bool is_correct_face(Faces_Ref faces) const;

    /// @brief На какой стороне расположена грань side для соседа, на
    /// структурированной сетке opposite_side(LEFT) = RIGHT и наоборот,
    /// а opposite_side(BOTTOM) = TOP и наоборот
    Side opposite_side(Side side) const;

    static_vector<Face_Ptr, 24> faces_list() const;

    static_vector<Face_Ptr, 24> faces() const;

    Faces_Ref faces(Side side) const;

    void set_faces(Side side, Faces_Ref faces);

    Vertex_Ref corner(Side side1, Side side2) const;

    Vertex_Ref left_top_vertex() const noexcept;

    Vertex_Ref left_bottom_vertex() const noexcept;

    Vertex_Ref right_bottom_vertex() const noexcept;

    Vertex_Ref right_top_vertex() const noexcept;

    /// @brief Возвращает все вершины (включая промежуточные на гранях),
    /// упорядоченные против часовой стрелки
    static_vector<Vertex *, 8> ordered_vertices() const;

    /// @brief Возвращает всех соседей
    static_vector<Cell*, 24> neighbors_raw() const;

    /// @brief Возвращает всех соседей
    static_vector<Cell::Ptr, 24> neighbors() const;

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Индексы.
     *
     *************************************************************************/

    /// @brief Проверяет, равны ли сеточные индексы и уровень данной ячейки и
    /// проверяемой.
    /// @param cell Ячейка для проверки идентичности.
    ///@return Вернёт true если x, y, lvl ячеек равны.
    bool equal_to(Cell::Ref cell) const noexcept;


    size_t level() const noexcept;

    void set_level(size_t lvl);

    size_t z_base() const noexcept;

    size_t z() const noexcept;

    void set_z(size_t z);


    IDType id() const noexcept;

    void set_id(IDType id);

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Иерархия ячеек.
     *
     *************************************************************************/

    Cell::Ptr parent() const noexcept;

    void set_parent(Cell::Ref parent);

    /// @brief Возвращает умный указатель на базовую ячейку, если ячейка
    /// сама является базовой, то возвращается её указатель.
    static Cell::Ptr origin(Cell::Ptr cell);


    const vector<Cell::Ptr> &children() const noexcept;

    /// @brief Дописывает листовых детей в вектор leaf_children, если дети отсутствуют,
    /// то дописывает указатель на себя
    static void leaf_children(Cell::Ref cell, vector<Cell::Ptr> &leaf_children);

    static_vector<Cell::Ptr, 4> children_by_side(Side side) const;

    /// @brief Возвращает дочернюю ячейку, которая является предковой для целевой
    /// @param nx -- временный костыль, число базовых ячеек по оси X
    Cell::Ref child_directed_to(const shared_ptr<Cell> &cell) const;

    void set_children(const vector<Cell::Ptr>& children);

    Cell::Ref get_child_by_local_id(size_t loc_z) const;

    static void set_child_by_local_id(Cell::Ref cell, size_t loc_z, Cell::Ref child);

    bool is_childless() const noexcept;

    void remove_children();


    static_vector<Cell *, 7> siblings() const;

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Данные.
     *
     *************************************************************************/

    int owner() const noexcept;

    void set_owner(int rank);

    bool is_local() const noexcept;

    bool is_remote() const noexcept;

    Data_Ref data_holder() const;

    double * data() const;

    void set_data(Data_Ref data);

    void copy_data_from(const Cell::Ptr &cell);


    AdaptationFlag adaptation_flag() const noexcept;

    void set_adaptation_flag(AdaptationFlag flag);

    bool can_coarse() const;

    double workload() const noexcept;

    void reset_workload();

    void add_workload(double workload);

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Геометрия.
     *
     * @brief Геометрические функции ячейки: объем, границы и т.д.
     *
     *************************************************************************/

    bool is_2D() const noexcept;

    bool is_3D() const noexcept;

    const Vector3d &center() const noexcept;

    const Vector3d &mass_center() const noexcept;

    Vector2d center_2d() const noexcept;

    double center_radius() const;

    double center_angle() const;

    double volume() const;

    double volume_as() const;

    /// @brief Вычисляет объемный интеграл Int (vec x - vec x0) dV
    Vector3d volume_as(const Vector3d & x0) const;

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Отладка
     *
     *************************************************************************/

    void print_dbg_info(bool neighbors_out = true) const;

    bool interesting() const;

    /** @} */

private:
    void initialize();


    /// @brief Грани -- основа ячейки, без граней ячейки быть не может.
    /// Любой конструктор ячейки требует грани.
    vector<Faces_Ptr> m_faces;

    /// @brief Центр ячейки
    Vector3d m_center;

    /// @brief Центр масс ячейки
    Vector3d m_mass_center;

    /// @brief Объем ячейки
    double m_volume;

    /// @brief Правильная ли ориентация граней, ориентация правильная, если
    /// внешняя нормаль к грани является внешней для данной ячейки
    vector<bool> m_or;

    /// @brief Уровень ячейки начиная с Cell::root_level().
    size_t m_level;

    /// @brief Положение ячейки на кривой, заполняющей пространство, если сетка
    /// имеет уровень адаптации до m_level
    /// Короче, в начале m_z выставляется некоторым образом в сеточном генераторе,
    /// для детей используются формулы z(lvl+1) = 4 z(lvl) + z_child_loc
    size_t m_z;

    /// @brief Уникальный идентификатор ячейке в иерархии.
    IDType m_id;

    /// @brief Указатель на родительскую ячейку.
    Cell::WPtr m_parent;

    /// @brief Указатели на дочерние ячейки (порядок важен).
    vector<Cell::Ptr> m_children;

    /// @{
    /// @name Данные ячейки.
    /// @brief Последующие поля являются опциональными, ячейка как
    /// геометрический/топологический элемент может существовать и без них.

    /// @brief Ранг процесса, который владеет ячейкой, для m_cells совпадает с
    /// mpi::rank(), для фейковых равен -1.
    int m_owner;

    ///< @brief Флаг адаптации.
    AdaptationFlag m_adaptation_flag;

    /// @brief Нагрузка на ячейку.
    double m_workload;

    /// @brief Указатель на данные ячейки
    Data_Ptr m_data;

    /// @} // Данные ячейки
};

/// @brief Структура ячейки для сериализации
struct ClearCell {
    IndexedObject::index_type index;
    IndexedObject::index_type faces[8];
    size_t id, z, level;
    int owner;
    bool fake;

    ClearCell() = default;

    ClearCell(Cell::Ref cell, bool is_fake);

    explicit ClearCell(const string& line);
};

/// @brief Вывод ячейки в ASCII поток
ostream& operator<<(ostream& out, const ClearCell& cell);
