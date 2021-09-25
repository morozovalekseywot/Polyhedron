//
// Created by 159-egi on 8/24/17.
//

#include "utils/memory/list_position.h"
#include "utils/memory/node_list.h"
#include "utils/memory/node.h"

ListPosition::ListPosition(NodeList *node_list, NodePtr node) :
        m_node(move(node)), m_next(nullptr), m_prev(nullptr) {

    if (!node_list) {
        throw runtime_error("Error in ListPosition constructor: nullptr node_list.");
    }
    if (!m_node) {
        throw runtime_error("Error in ListPosition constructor: nullptr node.");
    }
    m_node->m_positions[node_list] = this;
}

ListPosition::~ListPosition() {
    if (m_next)
        m_next->m_prev = m_prev;
    if (m_prev)
        m_prev->m_next = m_next;
}

void ListPosition::insert_before(ListPosition *prev) noexcept {
    if (m_prev) {
        m_prev->m_next = prev;
        prev->m_prev = m_prev;
    }
    m_prev = prev;
    prev->m_next = this;
}

ListPosition* ListPosition::insert_before(NodeList *node_list, NodePtr prev) noexcept {
    auto lp = new ListPosition(node_list, move(prev));
    insert_before(lp);
    return lp;
}

void ListPosition::insert_after(ListPosition *next) noexcept {
    if (m_next) {
        m_next->m_prev = next;
        next->m_next = m_next;
    }
    m_next = next;
    next->m_prev = this;
}

ListPosition* ListPosition::insert_after(NodeList *node_list, NodePtr next) noexcept {
    auto lp = new ListPosition(node_list, move(next));
    insert_after(lp);
    return lp;
}
