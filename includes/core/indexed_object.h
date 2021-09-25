//
// Created by 159-mrv on 05/29/18.
//

#ifndef INFINIMESH_CORE_INDEXED_OBJECT_H
#define INFINIMESH_CORE_INDEXED_OBJECT_H

/**
 * Макрос используется при пересылке индексов примитивов,
 * должен быть согласован с IndexedObject::index_type
 */
#define MPI_PRIMITIVE_INDEX MPI_SIZE_T


/**
 * @brief Абстрактный индекс.
 * Позволяет нумеровать экземпляры классов, которые его наследуют.
 */
class IndexedObject {
public:
    using index_type = size_t;

    static const index_type undefined_index = static_cast<index_type>(-1);

    IndexedObject() : m_index(undefined_index) { }

    index_type get_index() const { return m_index; }

    void set_index(index_type id) {  m_index = id; }

    void reset_index() { m_index = undefined_index; }

    bool is_index_undefined() const { return m_index == undefined_index; }

private:
    index_type m_index;
};


/**
 * Компоратор для проиндексированных объектов.
 * Позволяет сортировать массивы примитивов, унаследованных от IndexedObject,
 * а также хранить элементы в контейнерах типа set, map.
 */
struct CompareByIndex {

    using object_ptr = const shared_ptr<IndexedObject> &;
    using index = const IndexedObject::index_type &;

    bool operator()(object_ptr obj1, object_ptr obj2) const noexcept {
        return obj1->get_index() < obj2->get_index();
    }

    bool operator()(object_ptr obj, index ind) const noexcept {
        return obj->get_index() < ind;
    }

    bool operator()(index ind, object_ptr obj) const noexcept {
        return ind < obj->get_index();
    }
};


#endif //INFINIMESH_CORE_INDEXED_OBJECT_H
