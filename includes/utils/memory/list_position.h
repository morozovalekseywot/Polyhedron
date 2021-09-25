//
// Created by 159-egi on 8/23/17.
//

#ifndef AMR_LIST_POSITION_H
#define AMR_LIST_POSITION_H

#include <allstd.h>

class Node;
class ListPosition;
class NodeList;


class ListPosition {
public:
    using NodePtr = shared_ptr<Node>;

    /** Node, NodeList, ListPosition очень дружны */
    friend class Node;
    friend class NodeList;


    /** ListPosition содержит единственный Node, Node содержит много ListPosition */
    NodePtr m_node;

    /** Указатель на следующий ListPosition */
    ListPosition* m_next;

    /** Указатель на предыдущий ListPosition */
    ListPosition* m_prev;


    /**
     * @brief Конструктор класса
     * @param list Список в котором необходимо создать позицию
     * @param node Узел, на который ссылается позиция
     * @details Параметр node связывается с созданной ListPosition,
     * нет необходимости вызывать node->add_position(...) ещё раз.
     */
    ListPosition(NodeList* list, NodePtr node);

    /**
     * @brief Деструктор класса.
     * Удаляет ListPosition из цепочки, связывает все указатели
     */
    ~ListPosition();


    /** @defgroup list_position_insert_functions Функции вставки.
     *
     * @brief Функции в этой группе не поддерживают сохранность списка node_list.
     * То есть ListPositions вообще не знают,
     * что пренадлежат какому-либо списку.
     * Таким образом, не стоит вообще вызывать эти функции,
     * если вы не внутри node_list.
     *
     * @{
     */

    /**
     * @brief Вставляет Node после, связывает все указатели
     * @param list Список для вставки
     * @param prev Новый Node
     * @return Возвращает указатель на позицию,
     * созданную для нового узла prev
     */
    ListPosition* insert_before(NodeList* list, NodePtr prev) noexcept;

    /**
     * @brief Вставляет Node перед текущим, связывает все указатели
     * @param list Список для вставки
     * @param next Новый Node
     * @return Возвращает указатель на позицию,
     * созданную для нового узла next
     */
    ListPosition* insert_after(NodeList* list, NodePtr next) noexcept;

    /**
     * @brief Вставляет ListPosition после текущего, связывает все указатели
     * @param prev Новый ListPosition
     */
    void insert_before(ListPosition* prev) noexcept;

    /**
     * @brief Вставляет ListPosition перед текущим, связывает все указатели
     * @param next Новый ListPosition
     */
    void insert_after(ListPosition* next) noexcept;

    /** @} */ // list_position_insert

};

#endif //AMR_LIST_POSITION_H
