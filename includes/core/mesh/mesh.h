#pragma once

#include <allstd.h>

#include <core/indexed_object.h>
#include <core/cell/side.h>
#include <core/mesh/balancing_options.h>
#include <utils/math.h>
#include <utils/memory/node_list.h>

class Vertex;
class Face;
class FacesList;
class Cell;

enum class AdaptationFlag : int;
enum class Side : int;
enum class FaceFlag : int;

class BorderControl;
class GhostControl;

class Problem;
class Decomposition;
class GeometryGenerator;
class Configuration;


/**************************************************************************//**
 * @brief Класс создаёт и обслуживает сетку
 *
 * Данный класс предназначен для корректного создания и связывания
 * ячеек базовой сетки, гостов и граничных ячеек, а также граничных
 * условий. Также данный класс берёт на себя ответственность за
 * согласование и проведение адаптации сетки, балансировку нагрузки,
 * обмен данными с соседними процессами (включая их определение) и т.д.
 *****************************************************************************/
class Mesh {

    /// @brief Сокращенные именования внутри класса
    using Vertex_Ptr = shared_ptr<Vertex>;
    using Vertex_Ref = const shared_ptr<Vertex> &;
    using Face_Ptr   = shared_ptr<Face>;
    using Face_Ref   = const shared_ptr<Face> &;
    using Faces_Ptr  = shared_ptr<FacesList>;
    using Faces_Ref  = const shared_ptr<FacesList> &;
    using Cell_Ptr   = shared_ptr<Cell>;
    using Cell_Ref   = const shared_ptr<Cell> &;

    using GeometryGenerator_UPtr = unique_ptr<GeometryGenerator>;
    using Decomposition_UPtr = unique_ptr<Decomposition>;

    using BorderControl_UPtr = unique_ptr<BorderControl>;
    using GhostControl_UPtr  = unique_ptr<GhostControl>;


    /// @brief Геометрические типы адаптации
    enum class AdaptationCriterion : int {NONE, ALL, CASCADE };

public:
    using VertexVector = vector<Vertex_Ptr>;
    using FaceVector   = vector<Face_Ptr>;
    using CellVector   = vector<Cell_Ptr>;

    using UPtr = unique_ptr<Mesh>;


    /** @{ ****************************************************************//**
     *
     * @name Конструкторы.
     *
     *************************************************************************/

    /// @brief Осуществляет инициализацию базовых полей класса с
    /// использованием конфигурации. Затем вызывается функция build для
    /// непосредственной генерации сетки.
    /// @param config Экземпляр конфигурации
    explicit Mesh(const Configuration &config);

    /// @brief Возвращает уникальный указатель на сетку
    static Mesh::UPtr create(const Configuration &config);

    /// @brief Осуществляет инициализацию базовых полей класса с
    /// использованием конфигурации. Затем вызывается функция
    /// restart_from_checkpoint для непостредственной генерации сетки.
    /// @param config Экземпляр конфигурации
    /// @param checkpoint Путь к файлу с сохранением
    /// @warning Аргументы config и checkpoint должны быть согласованы,
    /// то есть сохранение checkpoint должно быть сделано с аналогичной
    /// конфигурацией
    Mesh(const Configuration &config, const string& checkpoint);

    /// @brief Возвращает уникальный указатель на сетку
    static Mesh::UPtr create(const Configuration &config, const string& checkpoint);

    /// @brief Деструктор
    ~Mesh();

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Методы доступа.
     *
     *************************************************************************/

    /// @brief Максимальный допустимый уровень адаптации
    size_t max_level();

    /// @brief Cписок всех ячеек локальной сетки
    NodeList::Ref cells();

    /// @brief Количество областей в разбиении по тредам
    uint n_chunks() const;

    /// @brief Полный спиок, альтернатива *Mesh::cells()
    NodeList::Part all_chunks();

    /// @brief Возвращает генератор геометрии
    const GeometryGenerator * geometry() const;

    /// @brief Вычисляет некоторую функцию на сетке
    /// @param func Некоторая функция от части ячеек сетки
    /// @return Время выполнения в секундах
    double map(const std::function<void(const NodeList::Part &)>& func);

    /// @brief Устанавливает размер DataHolder у всех ячеек равный cell_data_size
    void extend_cell_data(size_t cell_data_size);

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Основной функционал
     *
     *************************************************************************/

    /// @brief Выставляет флаги адаптации в ячейках, основываясь на алгоритме
    /// решателя. На соседних ячейках автоматически поддерживается баланс 2:1,
    /// предпочтение отдается ячейкам, которые хотят поделиться.
    /// @param problem Решатель, который отвечает за адаптацию и разделяет
    /// данные родительской ячейки между дочерними
    void adaptation(Problem* problem);

    /// @brief Синхронизация данных ячеек
    void sync_data();

    /// @brief Глобальная (сплошная) индексация листовых ячеек на всех процессах
    /// @param idx_begin Начальный индекс для ячеек данного процесса
    /// @param idx_end Конечный индекс для ячеек данного процесса (не достигается,
    /// n_cells_local = idx_end - idx_begin)
    /// @param n_cells_global Полное число ячеек на всех процессах
    void global_indexing(size_t& idx_begin, size_t& idx_end, size_t& n_cells_global);

    /// @brief Функция балансировки нагрузки
    /// @param proc_workload Нагрузка на текущий процесс, если не указана,
    /// то в качестве нагрузки используется число ячеек подсетки
    void load_balance(double proc_workload = 0.0);

    /// @brief Создание чекпоинта
    /// @param path Префикс имени файла секпоинта
    /// @param time Время
    /// @param binary Использовать бинарный формат
    void checkpoint(string path, double time, bool binary = true);

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Отладка
     *
     *************************************************************************/

    /// @brief Сохраняет все списки ячеек (m_cells, m_fake, m_ghost, m_border)
    /// в папку debug в формате VTU
    void dbg_write(Problem* problem, const string& prefix = "");

    /// @brief Выводит некоторую информацию о используемой сетке
    void print_info(const std::string& tab);

    /// @return Hash сетки - сумма хэшей данный каждой ячейки
    size_t hash();

    /** @} */

private:

    /** @{ ****************************************************************//**
     *
     * @name Инициализация, связывание ячеек
     *
     *************************************************************************/

    /// @brief Функция инициализирует ряд переменных, в частности
    /// определяются параметры адаптации, балансировки и сохранения файлов.
    /// Также в функции инициализируются такие переменные как m_geometry
    /// и m_decomposition, отвечающие за геометрию и деуомпозицию области,
    /// соответственно
    void read_configuration(const Configuration &config);

    /// @brief Последовательно вызывает все необходимые для создания сетки
    /// функции, начиная от создания вершин, заканчивая созданием гост-слоёв
    /// и связыванием ячеек сетки.
    void build(const Configuration &config);

    /// @brief Производит инициализацию ghost, border и fake ячеек.
    /// Для идентификации ячеек различных типов используется декомпозиция
    /// m_decomposition
    void initialize_other_cells(array<FaceFlag, 6> face_flags);

    /// @brief Соединяет ячейки по заданной стороне, подразумевается, что по
    /// заданной стороне существует ровно одна грань, в противном случае
    /// поведение неопределено, т.к. соединение происходит выборкой последней
    /// грани из массива граней данной ячейки по данному направлению.
    ///
    /// @param side Сторона, по которой необходимо связать ячейки
    /// @param cell Данная ячейка, по стороне которой происходит связь
    /// @param neighbor Соседняя ячейка, по оппозитной стороне которой происходит связь
    static void connect_cells(Side side, Cell_Ref cell, Cell_Ref neighbor);

    /// @brief Создает соседа для ячейки со стороны side
    /// @warning Данный метод подходит исключительно для создания фейковых
    /// ячеек!!! Метод не гарантирует создание корректного гост-слоя, а тем
    /// более базовой ячейки!!!
    Cell_Ptr create_neighbor(Cell_Ref cell, Side side);

    /// @brief Выполняет перестройку гост и бордер, которые могли измениться
    /// в ходе адаптации, балансировки или в другом случае
    void rebuild_ghosts();

    /// @brief Полностью удаляет связи ячеек внутри ghost-слоев, это необходимо
    /// для вменяемых операций дробления и огрубления ghost ячеек, также удаляет
    /// двойные стороны у ячеек в ghost-слое. Обязательно к применению после
    /// балансировки, когда может остаться иерархия
    void remove_ghost_connections();

    /// @brief Перемешивает список листовых ячеек и разбивает его на чанки
    void shuffle_and_split();

    /// @brief Z-индекс первой (нулевой) ячейки уровня выше level. При level = 0
    /// возвращает число базовых ячеек, при level = 1 возвращает сумму числа
    /// базовых ячеек и числа всех (возможных) ячеек первого уровня, и так далее
    size_t z_curve_offset_2d(size_t level) const;

    size_t z_curve_offset_3d(size_t level) const;

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Адаптация
     *
     *************************************************************************/

    /// @brief Адаптация базовой сетки исходя из геометрических соображений
    /// Вызывается один раз при создании сетки
    void geometric_adaptation(AdaptationCriterion criterion);

    /// @brief Общий критерий адатации, использует gtype и stype,
    /// а также учитывает ограничения на уровень
    /// @return Флаг адаптации
    AdaptationFlag adaptation_criterion(Problem* problem, Cell_Ref cell);

    /// @brief Устанавливает все флаги текущих ячеек и гостовых ячеек
    /// в значения AdaptationFlag::COARSE, если уровень позволяет (> 0)
    void reset_flags();

    /// @brief Устанавливает флаги в соответствии с критерием солвера
    /// Автоматически просиходит балансировка флагов в рамках данной части
    /// сетки
    void set_flags(Problem* problem);

    /// @brief Диагностическая функция, проверяет согласованность флагов,
    /// и возможность разбиения/огрубления
    void check_flags();

    /// @brief Установить ячейке флаг SPLIT и передать соседним ячейкам,
    /// чтобы они не отставали по уровням
    void set_need_to_split(Cell_Ref cell);

    /// @brief Установить ячейке флаг SPLIT и передать соседним ячейкам,
    /// чтобы они не отставали по уровням
    void set_need_to_split(Cell* cell);

    /// @brief Установить ячейке флаг NONE и передать соседним ячейкам,
    /// чтобы они не отставали по уровням
    void reset_adaptation_flag(Cell_Ref cell);

    /// @brief Установить ячейке флаг NONE и передать соседним ячейкам,
    /// чтобы они не отставали по уровням
    void reset_adaptation_flag(Cell* cell);

    /// @brief Установить ячейке флаг COARSE
    void set_need_to_coarse(Cell_Ref cell);

    /// @brief Функция вызывается после синхронизации флагов адаптации между
    /// процессами. Выполняет балансировку флагов данной части сетки, сохраняя
    /// флаги на гост-слоях неизменными.
    void balance_flags();

    /// @brief Применяет флаги адаптации
    /// @param problem Решатель, который отвечает за адаптацию
    void apply_flags(Problem* problem, bool debug = false);

    /// @param p_cell Ячейка, которую необходимо разбить
    /// @param ignore_neighbors Игнорировать соседей, если true, тогда связи с
    /// соседними ячейками не формируются, такой тип разбиения используется при
    /// восстановлении иерархии
    /// @param problem Решатель необходим для "умного" переноса данных от
    /// родительской ячейки к дочерним, при отсутствии решателя, данные
    /// полностью копируются в дочерние ячейки
    void split_cell(Cell_Ref p_cell, bool ignore_neighbors, Problem* problem = nullptr);

    /// @brief Соединяет дочерние ячейки с соседями cell, функция вызывается
    /// внутри split_cell после создания детей. На момент вызова функции дети
    /// связаны друг с другом, но ещё не имеют внешних связей с соседями cell
    /// @param cell Ячейка, которая подверглась разбиению
    void connect_children_with_neighbors(Cell_Ref cell);

    /// @brief Создает для двух дочерних ячеек по стороне side фейковых соседей
    /// и соединяет их друг с другом
    /// @param cell Ячейка, которая подверглась разбиению
    /// @param side Сторона ячейки у границы области
    void connect_children_with_fake_cells(Cell_Ref cell, Side side);

    /// @brief Соединяет две дочерних ячейки по стороне side с единственным
    /// соседом
    /// @param cell Ячейка, которая подверглась разбиению
    /// @param side Сторона ячейки у границы области
    void connect_children_with_single_cell(Cell_Ref cell, Side side);

    /// @brief Соединяет две дочерних ячейки по стороне side с двумя соседями
    /// с противоположной стороны
    /// @param cell Ячейка, которая подверглась разбиению
    /// @param side Сторона ячейки у границы области
    void connect_children_with_two_cells(Cell_Ref cell, Side side);

    /// @brief Вынимает родителя из всех списков, а в начало списков вставляет
    /// детей (выборочно, в зависимости от роли)
    /// @param cell Родительская ячейка
    void children_to_lists(Cell_Ref cell);

    /// @param cell Ячейка для объединения (должна иметь от 2 до 4 детей)
    /// @param problem Решатель необходим для объединения данных
    void coarse_cell(Cell_Ref cell, Problem* problem);

    /// @brief Вычисляет необходимый уровень адаптации по углу для радиуса r,
    /// при использовании каскадной адаптации, уровень ограничен переменной
    /// m_max_level.
    /// @param r Радиус
    /// @return Необходимый уровень адаптации
    size_t level_by_radius(double r);

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Балансировка
     *
     *************************************************************************/

    /// @brief Выполняет разделение сетки по процессам, в соответстии с
    /// разбиением выставляет поле Cell::m_owner для ячеек в списках m_cells,
    /// m_fake и m_ghosts
    /// @details Временная сложность: O(N) + O(splitter), где N - число
    /// ячеек подсетки.
    void compute_next_owners(double workload);

    vector<std::set<Cell_Ptr, CompareByIndex>> compute_cells_to_send() const;

    /// @brief Восстанавливает иерархию массиву ячеек, добавляя их корни в
    /// массив корневых ячеек.
    /// @param cells Массив ячеек для которых необходимо воссоздать иерархию
    void restore_hierarchy(CellVector& cells);

    /// @brief Номер гост-слоя в списке гостов, который соответствует данному
    /// рангу
    size_t add_neighbor_if_not_exist(int rank);

    /// @brief Восстанавливает иерархию от предковой ячейки до листовой
    void insert_child(Cell_Ptr parent, Cell_Ref child);

    /// @brief Данная функция необходимо во время балансировки нагрузки со
    /// стороны принимающего процесса. Проходит по списку ячеек и проверяет на
    /// равенство его индекса (z,lvl), если такой ячейки не существует
    /// в списке - добавляет её в этот список.
    /// @param cell Ячейка к добавлению
    /// @param cells Список, в котором будет произведён поиск индекса ячейки
    void push_cell_if_not_equal(Cell_Ref cell, NodeList::Ptr cells);

    /// @brief Содержит ли данная часть сетки базовую ячейку для заданой
    bool base_cell_exist_for(Cell_Ref cell);

    /// @brief Добавить базовую ячейку
    void insert_base_cell(Cell_Ref cell);

    /// @brief Взять базовую ячейку по индексу
    Cell_Ptr base_cell(size_t z);

    /// В случае, когда после балансировки ячейка из основного списка переходит
    /// в гост-слой, у соответственной базовой ячейки сохраняется вся иерархия.
    /// Это может помешать огрублению ячеек гост-слоя, поскольку базовая ячейка
    /// будет содержать множество "неактуальных" детей (не включенных в m_cells,
    /// m_ghosts или m_fake). Таким образом, гост-слой необходимо поддерживать
    /// на минимальном уровне адаптации, чем и занимается данная функция.
    /// @brief Выполняет очистку иерархии гост-слоя после балансировки.
    void clean_ghosts();

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Рестарт
     *
     *************************************************************************/

    /// @brief Вспомогательная функция, для согласования версий чекпоинтов
    double checkpoint_converter(string full_path);

    /// @brief Восстановление из чекпоинта. Заменяет функцию build() внутри
    /// конструктора, если программа запущена с чекпоинта
    double restart_from_checkpoint(string full_path);

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Коммуникация
     *
     *************************************************************************/

    /// @brief Функция вызывается для страховки, высылает количество ячеек в
    /// граничных слоях, дожидается приёма количества ячеек в гост-слоях. Если
    /// количество ячеек в нашем гост-слое не совпадёт в количеством ячеек в
    /// удалённом граничном слое - бросит исключение.
    void sync_cells_count();

    /// @brief Осуществляет обмен флагами адаптации до тех пор, пока
    /// функция ghost_flags_changed на всех процессах не начнёт возвращать
    /// false. После синхронизации флагов адаптации приграничных ячеек
    /// происходит вызов функции балансировки флагов balance_flags.
    void sync_adaptation();

    /// @brief Синхроинизирует флаги адаптации
    void sync_flags();

    /// @brief Проходит по всем гост-слоям и проверяет изменились ли у них
    /// флаги адаптации после последней синхронизации
    /// @return True, если флаги в гост-слое изменились после синхронизации
    bool ghost_flags_changed();

    /// @brief Сквозная (не сплошная) индексация всех вершин, сторон и ячеек,
    /// включает нумерацию базовых и фейковых ячеек,
    /// нумерация ghost ячеек синхронизуется.
    /// Функция используется при рестарте и балансировке.
    /// Аргументы являются выходными параметрами.
    /// @param vertices Вектор всех вершин сетки (m_cells + m_fake)
    /// @param faces    Вектор всех простых граней (m_cells + m_fake)
    /// @param cells    Вектор всех ячеек (m_cells + m_fake)
    void global_indexing(VertexVector& vertices, FaceVector& faces, CellVector& cells);

    /// @brief Сквозная (не сплошная) индексация всех вершин, сторон и ячеек,
    /// включает нумерацию базовых и фейковых ячеек,
    /// нумерация ghost ячеек синхронизуется.
    /// Функция используется при рестарте и балансировке.
    /// Аргументы являются выходными параметрами.
    /// @param vertices Контейнер содержит все вершины сетки (m_cells + m_fake)
    /// @param faces    Контейнер содержит все простые грани (m_cells + m_fake)
    /// @param cells    Контейнер содержит все ячейки (m_cells + m_fake)
    void global_indexing(std::map<IndexedObject::index_type, Vertex_Ptr> &vertices,
                         std::map<IndexedObject::index_type, Face_Ptr> &faces,
                         std::map<IndexedObject::index_type, Cell_Ptr> &cells);

    /// @brief Синхронизирует индексацию Ghost слоев.
    /// При выборе из двух возможных индексов (для вершин и граней)
    /// предпочтение отдается меньшему
    void sync_indices();

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Отладка
     *
     *************************************************************************/

    void dbg_write_seq();

    void dbg_check_duplicates();

    void dbg_write_connectivity(size_t step);

    void dbg_write_ghost_connectivity(size_t step);

    /** @} */


private:

    /// @brief Генератор геометрии сетки
    GeometryGenerator_UPtr m_geometry;

    /// @brief Класс, осуществляющий декомпозицию сетки
    Decomposition_UPtr m_decomposition;

    /**
     * Списки ячеек данной части сетки
     *
     * m_base_cells - локальные базовые ячейки
     * m_cells      - локальные листовые ячейки
     * m_fake       - локальные фековые ячейки (т.е. смежные с m_cells)
     * m_border     - вектор границ данной сетки
     * m_ghosts     - ячейки, смежные к m_border
     */
    unordered_map<size_t, Cell_Ptr> m_base_cells;

    NodeList::Ptr m_cells;
    NodeList::Ptr m_fake;

    vector<BorderControl_UPtr> m_border;
    vector< GhostControl_UPtr> m_ghosts;

    /// @brief Максимальный уровень адаптации
    size_t m_max_level;

    /// @brief Геометрический тип адаптации.
    AdaptationCriterion m_adaptation_type;

    /// @brief Ипользовать критрий адаптаций решателя
    bool m_solver_adaptation;

    /// @brief Опции балансировки.
    BalancingOptions m_balancing_options;

    /// @brief Параллельность по тредам
    uint m_n_chunks;
};
