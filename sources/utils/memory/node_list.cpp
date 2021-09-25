//
// Created by 159-egi on 8/24/17.
//

//#include "control/parallel_controller.h"
#include <random>
#include "core/cell/cell.h"

#include "utils/memory/node.h"
#include "utils/memory/list_position.h"
#include "utils/memory/node_list.h"

NodeList::NodeList(int rank, ListRole role)
  : m_head(nullptr),
    m_tail(nullptr),
    m_size(0),
    m_rank(rank),
    m_role(role),
    m_parts({m_head, nullptr}){
}

NodeList::~NodeList() {
    clear();
}

void NodeList::remove_parts() noexcept {
    m_parts.resize(2);
    m_parts[0] = m_head;
    m_parts[1] = nullptr;
}

void NodeList::push_back(Node_Ptr new_node) {
    if (contain(new_node.get())) {
        return;
    }

    m_mutex.lock();

    if (m_tail) {
        m_tail = m_tail->insert_after(this, move(new_node));
    }
    else {
        m_tail = new ListPosition(this, move(new_node));
        if (!m_head) {
            m_head = m_tail;
        }
    }
    remove_parts();

    m_mutex.unlock();

    ++m_size;
}

void NodeList::push_forward(Node_Ptr new_node) {
    if (contain(new_node.get())) {
        return;
    }

    m_mutex.lock();

    if (m_head) {
        m_head = m_head->insert_before(this, move(new_node));
    }
    else {
        m_head = new ListPosition(this, move(new_node));
        if (!m_tail) {
            m_tail = m_head;
        }
    }

    remove_parts();

    m_mutex.unlock();

    ++m_size;
}

bool NodeList::contain(Node_Ptr node) const {
    if (!node)
        return false;
    auto list_ptr = const_cast<NodeList*>(this);
    return node->m_positions.find(list_ptr) != node->m_positions.end();
}

bool NodeList::contain(Node* node) const {
    if (!node)
        return false;
    auto list_ptr = const_cast<NodeList*>(this);
    return node->m_positions.find(list_ptr) != node->m_positions.end();
}

void NodeList::erase(Node_Ptr node) {
    erase(node.get());
}

void NodeList::erase(Node* node) {
    m_mutex.lock();

    if (node == m_head->m_node.get()) {
        m_head = m_head->m_next;
    }
    if (node == m_tail->m_node.get()) {
        m_tail = m_tail->m_prev;
    }

    auto it = node->m_positions.find(this);

    if (it != node->m_positions.end()) {
        delete it->second;
        node->m_positions.erase(it);
    }
    else {
        throw runtime_error("Error! Node to delete not found");
    }

    remove_parts();

    m_mutex.unlock();

    --m_size;
}

void NodeList::clear() {
    m_mutex.lock();
    auto it = m_head;
    while (it != nullptr) {
        auto lp_to_delete = it;
        it = it->m_next;
        lp_to_delete->m_node->m_positions.erase(this);
        delete lp_to_delete;
    }
    m_tail = nullptr;
    m_head = nullptr;

    remove_parts();

    m_mutex.unlock();

    m_size = 0;
}

bool NodeList::empty() const noexcept {
    return m_head == nullptr;
}

void NodeList::set_role(ListRole role) noexcept {
    m_role = role;
}

ListRole NodeList::role() const noexcept {
    return m_role;
}

size_t NodeList::size() const noexcept {
    return m_size;
}

void NodeList::shuffle(bool true_random) {
    // Два статических генератора
    static std::mt19937_64 true_generator(time(0));
    static std::mt19937_64 false_generator(0);

    // Выбираем необходимый
    std::mt19937_64& generator = true_random ? true_generator : false_generator;

    // Равномерное распределение
    std::uniform_int_distribution<int> uniform(0, std::numeric_limits<int>::max());

    vector<ListPosition*> elements;
    elements.reserve(m_size);

    ListPosition* lp = m_head;
    while (lp != nullptr) {
        elements.emplace_back(lp);
        lp = lp->m_next;
    }

    // Khnut shuffle
    for(auto i = int(m_size - 1); i > 0; --i) {
        int j = uniform(generator) % (i + 1);

        auto temp = elements[i];
        elements[i] = elements[j];
        elements[j] = temp;
    }

    m_head = elements[0];
    m_head->m_prev = nullptr;
    m_head->m_next = m_size > 1 ? elements[1] : nullptr;

    m_tail = elements[m_size - 1];
    m_tail->m_next = nullptr;
    m_tail->m_prev = m_size > 1 ? elements[m_size - 2] : nullptr;

    for(size_t i = 1; i < m_size - 1; ++i) {
        elements[i]->m_prev = elements[i - 1];
        elements[i]->m_next = elements[i + 1];
    }

    remove_parts();
}

int NodeList::rank() const noexcept {
    return m_rank;
}


// ============================================================================
// Простой итератор
// ============================================================================

NodeList::iterator::iterator(ListPosition* pos) noexcept
        : m_position(pos) {

}

void NodeList::iterator::operator++() noexcept {
    m_position = m_position->m_next;
}

void NodeList::iterator::operator--() noexcept {
    m_position = m_position->m_prev;
}

Cell::Ptr NodeList::iterator::operator*() const {
    return std::static_pointer_cast<Cell>(m_position->m_node);
}

Cell::Ptr NodeList::iterator::operator->() const {
    return std::static_pointer_cast<Cell>(m_position->m_node);
}

bool NodeList::iterator::operator!=(const NodeList::iterator& it) const noexcept {
    return m_position != it.m_position;
}

NodeList::iterator NodeList::begin() const noexcept {
    return iterator(m_head);
}

NodeList::iterator NodeList::end() const noexcept {
    return iterator(nullptr);
}

void NodeList::split(uint nparts) {
    // Размеры маленького и большого чанков
    size_t lit_chunk = m_size / nparts;
    size_t big_chunk = m_size / nparts + 1;

    // Количество маленьких и больших чанков
    size_t lit_count = big_chunk * nparts - m_size;
    size_t big_count = nparts - lit_count;

    m_parts.clear();
    auto position = m_head;
    for (size_t i = 0; i < big_count; ++i) {
        m_parts.push_back(position);
        for (size_t j = 0; j < big_chunk; ++j) {
            position = position->m_next;
        }
    }
    for (size_t i = 0; i < lit_count; ++i) {
        m_parts.push_back(position);
        for (size_t j = 0; j < lit_chunk; ++j) {
            position = position->m_next;
        }
    }
    m_parts.push_back(position);
}

uint NodeList::n_parts() const noexcept {
    return static_cast<uint>(m_parts.size() - 1);
}

NodeList::Part NodeList::part(uint idx) const noexcept {
    return { m_parts[idx], m_parts[idx + 1], idx };
}

NodeList::Part NodeList::all_parts() const noexcept {
    return { m_head, nullptr, 0 };
}

NodeList::Part::Part(ListPosition *begin, ListPosition *end, uint id) noexcept
    : m_begin(begin), m_end(end), m_id(id) {
}

NodeList::iterator NodeList::Part::begin() const noexcept {
    return NodeList::iterator(m_begin);
}

NodeList::iterator NodeList::Part::end() const noexcept {
    return NodeList::iterator(m_end);
}

uint NodeList::Part::part_id() const noexcept {
    return m_id;
}