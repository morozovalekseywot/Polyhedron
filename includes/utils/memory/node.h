//
// Created by 159-egi on 8/23/17.
//

#ifndef AMR_NODE_H
#define AMR_NODE_H

#include <allstd.h>

class ListPosition;
enum class ListRole;
class NodeList;

/// @brief Элемент NodeList, может содержаться одновременно в нескольких NodeList.
/// @details Абстрактный класс, любой наследник может храниться в NodeList
/// (на данный момент его наследует только Cell)
class Node {

    /** Node, NodeList, ListPosition очень дружны */
    friend NodeList;
    friend ListPosition;

public:

    /// @brief Конструктор по умолчанию
    Node() = default;

    /// @brief Виртуральный деструктор
    virtual ~Node() = default;


    /// @defgroup node_simple_get Простые get методы.
    /// @{

    /// @brief Имеет ли узел роль role
    bool has_role(ListRole role) const noexcept;

    /// @brief Список NodeList, которому принадлежит узел.
    /// @throw runtime_error в случае, если узел принадлежит не одному списку.
    NodeList* get_list() const;

    /// @brief Множество списков, которым принадлежит узел
    vector<NodeList*> get_lists() const noexcept;

    /// @} // node_simple_get

    void remove_from_all_lists();

private:
    std::map<NodeList*, ListPosition*> m_positions;
};

#endif //AMR_NODE_H
