//
// Created by 159-mrv on 9/13/18.
//

#ifndef INFINIMESH_FACES_LIST_H
#define INFINIMESH_FACES_LIST_H

#include <allstd.h>

class Vertex;
class Face;
class Cell;

/**************************************************************************//**
 * @brief Класс для хранения сложных граней.
 *
 * @details Внутренние грани упорядочиваются таким образом, чтобы конец первой
 * грани совпадал с началом второй. Все операции с классом сохраняют выполнение
 * этого свойства.
 *****************************************************************************/
class FacesList : public std::enable_shared_from_this<FacesList> {

    /// @{ @brief Сокращения, используемые в описании класса.
    using Vertex_Ptr = shared_ptr<Vertex>;
    using Vertex_Ref = const shared_ptr<Vertex> &;
    using Face_Ptr   = shared_ptr<Face>;
    using Face_Ref   = const shared_ptr<Face> &;
    using Cell_Ptr   = shared_ptr<Cell>;
    using Cell_Ref   = const shared_ptr<Cell> &;
    /// @}


public:
    /// @brief Сокращенное именование shared_ptr на FacesList.
    using Ptr = shared_ptr<FacesList>;

    /// @brief Константная ссылка на shared_ptr для передачи в функции.
    using Ref = const shared_ptr<FacesList> &;



    /** @{ ****************************************************************//**
     *
     * @name Конструторы.
     *
     *************************************************************************/

    FacesList() = default;

    explicit FacesList(Face_Ref face);

    FacesList(Face_Ref face_1, Face_Ref face_2);

#if DIM3
    FacesList(Face_Ref face_1, Face_Ref face_2, Face_Ref face_3, Face_Ref face_4);
#endif

    static inline FacesList::Ptr create() {
        return make_shared<FacesList>();
    }

    static inline FacesList::Ptr create(Face_Ref face) {
        return make_shared<FacesList>(face);
    }

    static inline FacesList::Ptr create(Face_Ref face_1, Face_Ref face_2) {
        return make_shared<FacesList>(face_1, face_2);
    }

#if DIM3
    static inline FacesList::Ptr create(Face_Ref face_1, Face_Ref face_2,
                                        Face_Ref face_3, Face_Ref face_4) {
        return make_shared<FacesList>(face_1, face_2, face_3, face_4);
    }
#endif

    static inline FacesList::Ptr create(static_vector<Face_Ptr, 4> faces) {
        if (faces.size() < 3) {
            return make_shared<FacesList>(faces[0], faces[1]);
        }
#if DIM3
        else {
            return make_shared<FacesList>(faces[0], faces[1], faces[2], faces[3]);
        }
#endif
    }

    /** @} */



    /** @{ ****************************************************************//**
     *
     * @name Численные характеристики.
     *
     *************************************************************************/

    /// @brief Возвращает среднее арифметическое центров подграней
    Vector3d center() const noexcept;

    /// @brief Внешняя нормаль к сложной грани (среднее арифметическое)
    Vector3d outward_normal() const noexcept;

    /// @brief Площадь/длина сложной грани
    double area() const noexcept;

    /* @} */



    /** @{ ****************************************************************//**
     *
     * @name Структура грани.
     * Методы доступа к вершинам, подграням и смежным ячейкам.
     *
     *************************************************************************/

    /// @brief Возвращает shared_ptr на вершину по индексу так, как если бы
    /// данная грань была простой
    Vertex_Ref vertex(int idx) const;

    /// @brief Возвращает указатель на вершину по индексу так, как если бы
    /// данная грань была простой
    Vertex* vertex_raw(int idx) const;

    /// @brief Является ли грань одинарной?
    bool is_simple() const noexcept;

    /// @brief Является ли грань сложной?
    bool is_complex() const noexcept;

    /// @brief Возвращает число подграней.
    size_t size() const noexcept;

    /// @brief Возвращает ссылку на подгрань с номером idx.
    Face_Ref at(size_t idx) const noexcept;

    /// @brief Возвращает простую подгрань
    Face_Ref single()  const noexcept;

    /// @brief Добавляет всех соседей ячейки cell в вектор.
    /// @param cell Целевая ячейка, для которой ищутся соседи.
    /// @param neighbors Вектор ячеек, в который добавляются соседи.
    void append_neighbors(const Cell* cell, static_vector<Cell*, 24>& neighbors) const;

    /// @brief Добавляет всех соседей ячейки cell в вектор.
    /// @param cell Целевая ячейка, для которой ищутся соседи.
    /// @param neighbors Вектор ячеек, в который добавляются соседи.
    void append_neighbors(const Cell* cell, static_vector<Cell_Ptr, 24> & neighbors) const;

    /** @} */

private:
    vector<Face_Ptr> m_faces;
};

#endif //INFINIMESH_FACES_LIST_H
