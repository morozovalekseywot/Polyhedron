//
// Created by 159-egi on 8/23/17.
//

#ifndef AMR_LINKED_LIST_H
#define AMR_LINKED_LIST_H

#include <allstd.h>

class Cell;
class Node;
class ListPosition;

enum class ListRole {
    GHOST, BORDER, MESH
};

class NodeList {
    using Node_Ptr = shared_ptr<Node>;
    using Node_Ref = const shared_ptr<Node> &;

    /** Node, NodeList, ListPosition очень дружны */
    friend Node;
    friend ListPosition;

public:
    using Ptr = shared_ptr<NodeList>;
    using Ref = const shared_ptr<NodeList> &;

    explicit NodeList(int rank, ListRole role = ListRole::MESH);
    ~NodeList();

    /**
     * Копирование и перемещение запрещено,
     * Node используют указатель на список
     */
    NodeList(const NodeList& node_list) = delete;
    NodeList(const NodeList&& node_list) = delete;


    static inline NodeList::Ptr create(int rank, ListRole role = ListRole::MESH) {
        return make_shared<NodeList>(rank, role);
    }

    void push_back   (Node_Ptr new_node);
    void push_forward(Node_Ptr new_node);

    void erase(Node_Ptr node);
    void erase(Node* node);

    void clear();
    bool contain(Node_Ptr node) const;
    bool contain(Node* node) const;
    bool empty() const noexcept;
    size_t size() const noexcept;

    /// @brief Перемешивает список. Данная операция может быть использована
    /// для выравнивания нагрузки между тредами.
    /// @param true_random Использовать полность случайный рандомайзер, или
    /// засеивать генератор одним числом.
    void shuffle(bool true_random = false);

    int rank() const noexcept;
    ListRole role() const noexcept;
    void set_role(ListRole role) noexcept;


/** @brief Простой двунаправленный итератор по списку ячеек.
     * @example Примеры использования в циклах:
     * NodeList::Ptr cells;
     * ...
     * for(auto cell: *cells) {
     *     // cell имеет тип Cell::Ptr
     *     ...
     * }
     * for(auto it = cells->begin(); it != cells->end(); ++it) {
     *     Cell::Ptr cell = *it;
     *     double x = it->faces();
     *     ...
     * }     *
     * @warning Операции с изменением списка (insert, erase ...)
     * могут привести к неликвидности итератора
     */
    class iterator {
        friend NodeList;

    public:
        void operator++() noexcept;
        void operator--() noexcept;
        shared_ptr<Cell> operator*() const;
        shared_ptr<Cell> operator->() const;
        bool operator!=(const NodeList::iterator& it) const noexcept;

    public:
        explicit iterator(ListPosition* pos) noexcept;
        ListPosition* m_position;
    };
    iterator begin() const noexcept;
    iterator end() const noexcept;


    void split(uint nparts);

    class Part {
        friend NodeList;

    public:
        iterator begin() const noexcept;
        iterator end() const noexcept;
        uint part_id() const noexcept;

    private:
        Part(ListPosition* begin, ListPosition* end, uint id) noexcept;
        ListPosition* m_begin, *m_end;
        uint m_id;
    };

    uint n_parts() const noexcept;

    Part part(uint idx) const noexcept;

    Part all_parts() const noexcept;

private:

    void remove_parts() noexcept;

    /** Указатель на начало списка */
    ListPosition* m_head;

    /** Указатель на конец списка */
    ListPosition* m_tail;

    /** Размер списка */
    std::atomic_size_t m_size;

    /** Целевой ранк для передачи/приёма */
    int m_rank;

    /** Роль списка: MESH, BORDER или GHOST */
    ListRole m_role;

    /** Штуковина для потокобезопасности */
    std::mutex m_mutex;

    /** Разбиение списка на части, первые элементы частей + nullptr в качестве последнего элемента */
    vector<ListPosition*> m_parts;
};

#endif //AMR_LINKED_LIST_H
