//
// Created by 159-mrv on 9/13/18.
//

#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/cell.h>
#include <utils/geom.h>

FacesList::FacesList(Face::Ref face) {
    m_faces = {face};
}

FacesList::FacesList(Face::Ref face_1, Face::Ref face_2) {
    if (!face_1 || !face_2) {
        throw runtime_error("FacesList::FacesList(face1, face2) error: Attempt to create FacesList without faces");
    }

    // Направление нормалей должно совпадать
    if (face_1->outward_normal().dot(face_2->outward_normal()) < 0.0) {
        face_2->reverse();
    }

    // Правильный порядок
    if (face_1->vertex(1) == face_2->vertex(0)) {
        m_faces = {face_1, face_2};
        return;
    }
    // Обратный порядок
    if (face_2->vertex(1) == face_1->vertex(0)) {
        m_faces = {face_2, face_1};
        return;
    }

    // Простые грани не имеют общих вершин,
    // Возможно, необходимо склеить вершины
    double distance = Vertex::distance2(face_1->vertex(1), face_2->vertex(0));
    if (distance < 1.0e-6 * face_1->area()) {
        face_2->set_vertex(0, face_1->vertex(1));
        m_faces = {face_1, face_2};
        return;
    }
    distance = Vertex::distance2(face_2->vertex(1), face_1->vertex(0));
    if (distance < 1.0e-6 * face_1->area()) {
        face_2->set_vertex(1, face_1->vertex(0));
        m_faces = {face_2, face_1};
        return;
    }

    // Случай добавления эквивалентной стороны,
    // иногда возникает при рестарте, причину выяснять лень
    distance = Vertex::distance2(face_1->vertex(0), face_2->vertex(0)) +
               Vertex::distance2(face_1->vertex(1), face_2->vertex(1));
    if (distance < 1.0e-6 * face_1->area()) {
        face_2->set_vertex(0, face_1->vertex(0));
        face_2->set_vertex(1, face_1->vertex(1));
        m_faces = {face_1};
        return;
    }

    std::cerr << "    FacesList::create(face1, face2) error: Two simple faces have no shared vertex\n";
    std::cerr << "      face1: " << to_string(face_1) << "\n";
    std::cerr << "      face2: " << to_string(face_2) << "\n";

    throw runtime_error("FacesList::create(face1, face2) error: Two simple faces have no shared vertex");
}

#if DIM3
FacesList::FacesList(Face::Ref face_1, Face::Ref face_2, Face::Ref face_3, Face::Ref face_4) {
    if (!face_1 || !face_2 || !face_3 || !face_4) {
        throw runtime_error(
                "FacesList::FacesList(face1, face2, face3, face4) error: Attempt to create FacesList without faces");
    }

    array<Face::Ptr, 4> face = {face_1, face_2, face_3, face_4};
    array<Vector3d, 4> n;
    for (int i = 0; i < 4; ++i) {
        n[i] = face[i]->outward_normal();

        // Направление нормалей должно совпадать
        if (n[i].dot(n[0]) < 0.0) {
            face[i]->reverse();
            n[i] *= -1.0;
        }
    }

    // Найдем систему координат, связанную со сложной гранью
    // Начало координат - центр сложной грани
    // Ox - вектор в направлении от центра первой грани к общему
    // центру, таким образом полярный угол центра первой грани
    // в плоскости Oxy равен -Pi
    // Oz - внешняя нормаль сложной грани
    Vector3d center = Vector3d::Zero();
    Vector3d Ox = Vector3d::Zero();
    Vector3d Oy = Vector3d::Zero();
    Vector3d Oz = Vector3d::Zero();
    for (int i = 0; i < 4; ++i) {
        Oz += n[i];
        center += face[i]->center();
    }
    Oz.normalize();
    center /= 4.0;

    Ox = center - face[0]->center();
    Ox = Ox - Oz * Ox.dot(Oz);
    Ox.normalize();

    Oy = Oz.cross(Ox);

    // Углы центров граней в плоскости Oxy
    array<double, 4> phi = {
            -M_PI,
            atan2(Oy.dot(face[1]->center() - center), Ox.dot(face[1]->center() - center)),
            atan2(Oy.dot(face[2]->center() - center), Ox.dot(face[2]->center() - center)),
            atan2(Oy.dot(face[3]->center() - center), Ox.dot(face[3]->center() - center))
    };

    // Сортируем face_2, face_3, face_4 пузырьком по увеличению угла центра грани
    // в плоскости Oxy
    if (phi[1] > phi[2]) {
        std::swap(face[1], face[2]);
        std::swap(phi[1], phi[2]);
    }
    if (phi[2] > phi[3]) {
        std::swap(face[2], face[3]);
        std::swap(phi[2], phi[3]);
    }
    if (phi[1] > phi[2]) {
        std::swap(face[1], face[2]);
        std::swap(phi[1], phi[2]);
    }


    // Теперь гарантируется, что направление обхода граней face_1, face_2, face_3, face_4
    // совпадает с направлением охода вершин в любой из граней

    // Характерная длина
    double L = Vertex::distance2(face[0]->vertex(1), face[0]->vertex(0));

    // Вращаем грани до необходимого расположения
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            double d1 = Vertex::distance2(face[0]->vertex(i), face[2]->vertex(j));
            double d2 = Vertex::distance2(face[1]->vertex(i), face[3]->vertex(j));
            if (d1 < 1.0e-6 * L) {
                face[0]->rotate(i - 2);
                face[2]->rotate(j);
            }
            if (d2 < 1.0e-6 * L) {
                face[1]->rotate(i + 1);
                face[3]->rotate(j - 1);
            }
        }
    }

    // Склеиваем точки, если это необходимо
    Vertex::Ptr center_vertex = face[0]->vertex(2);

    for(int i = 0; i < 4; ++i) {
        int j = (i + 1) % 4;
        int k = (i + 2) % 4;
        if (Vertex::distance2(face[i]->vertex(j), face[j]->vertex(i)) > 1.0e-6 * L ||
            Vertex::distance2(face[i]->vertex(k), center_vertex) > 1.0e-6 * L) {
            std::cerr
                    << "    FacesList::create(face1, face2, face3, face4) error: Four simple faces have no shared vertex\n";
            std::cerr << "      face1: " << to_string(face[0]) << "\n";
            std::cerr << "      face2: " << to_string(face[1]) << "\n";
            std::cerr << "      face3: " << to_string(face[2]) << "\n";
            std::cerr << "      face4: " << to_string(face[3]) << "\n";

            throw runtime_error(
                    "FacesList::create(face1, face2, face3, face4) error: Four simple faces have no shared vertex");
        }

        face[j]->set_vertex(i, face[i]->vertex(j));
        face[i]->set_vertex(k, center_vertex);
    }

    m_faces = {face[0], face[1], face[2], face[3]};
}
#endif

Vertex::Ref FacesList::vertex(int idx) const {
    return is_simple() ? m_faces[0]->vertex(idx) : m_faces[idx]->vertex(idx);
}

Vertex * FacesList::vertex_raw(int idx) const {
    return is_simple() ? m_faces[0]->vertex_raw(idx) : m_faces[idx]->vertex_raw(idx);
}

Vector3d FacesList::center() const noexcept {
    Vector3d cntr = { 0.0, 0.0, 0.0 };
    for(Face::Ref face: m_faces) {
        cntr += face->center();
    }
    return cntr / m_faces.size();
}

Vector3d FacesList::outward_normal() const noexcept {
    Vector3d norm = { 0.0, 0.0, 0.0 };
    for(Face::Ref face: m_faces) {
        norm += face->outward_normal();
    }
    norm.normalize();
    return norm;
}

double FacesList::area() const noexcept {
    double res = 0.0;
    for(auto& f: m_faces) {
        res += f->area();
    }
    return res;
}

size_t FacesList::size() const noexcept {
    return m_faces.size();
}

bool FacesList::is_simple() const noexcept {
    return m_faces.size() == 1;
}

bool FacesList::is_complex() const noexcept {
    return m_faces.size() == 2;
}

Face::Ref FacesList::at(size_t idx) const noexcept {
    return m_faces[idx];
}

Face::Ref FacesList::single()  const noexcept {
    return m_faces.front();
}

void FacesList::append_neighbors(const Cell* cell, static_vector<Cell*, 24>& neighbors) const {
    for (Face::Ref face: m_faces) {
        if (face->neighbor_exist(cell)) {
            neighbors.emplace_back(face->neighbor_raw(cell));
        }
    }
}

void FacesList::append_neighbors(const Cell* cell, static_vector<Cell::Ptr, 24>& neighbors) const {
    for (Face::Ref face: m_faces) {
        if (face->neighbor_exist(cell)) {
            neighbors.emplace_back(face->neighbor(cell));
        }
    }
}