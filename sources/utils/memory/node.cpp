//
// Created by 159-egi on 8/23/17.
//

#include "utils/memory/list_position.h"
#include "utils/memory/node_list.h"
#include "utils/memory/node.h"

bool Node::has_role(ListRole role) const noexcept {
    for(const auto& position: m_positions) {
        if (position.first->role() == role) {
            return true;
        }
    }
    return false;
}

NodeList* Node::get_list() const {
    if (m_positions.size() == 1) {
        return m_positions.begin()->first;
    }
    else {
        throw runtime_error("ERROR! Attempt to get list from Node which has many lists.");
    }
}

vector<NodeList*> Node::get_lists() const noexcept {
    vector<NodeList*> list(m_positions.size());
    size_t i = 0;
    for(auto& kek: m_positions) {
        list[i++] = kek.first;
    }
    return list;
}


void Node::remove_from_all_lists() {
    for(auto it = m_positions.begin(); it != m_positions.end(); ) {
        auto list = it->first;

        if (this == list->m_head->m_node.get()) {
            list->m_head = list->m_head->m_next;
        }
        if (this == list->m_tail->m_node.get()) {
            list->m_tail = list->m_tail->m_prev;
        }

        delete it->second;
        it = m_positions.erase(it);

        list->remove_parts();
        --list->m_size;
    }
}