//
// Created by 159-mrv on 4/25/18.
//

#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/cell/cell.h>
#include <core/layer/border_control.h>
#include <core/layer/ghost_control.h>

#include <core/mesh/mesh.h>


// ===============================================================================
// =================== 03 - Часть отвечающая за коммуникацию =====================
// ===============================================================================

void Mesh::sync_cells_count() {
    // Синхронизация количества ячеек
    for (auto &border: m_border) { border->send_cells_count(); }
    for (auto &ghost:  m_ghosts) { ghost->recv_cells_count(); }

    for (auto &border: m_border) { border->wait_cells_count(); }
    for (auto &ghost:  m_ghosts) { ghost->wait_cells_count(); }
}

void Mesh::sync_data() {
    sync_cells_count();

    for (auto &border: m_border) { border->send_data(); }
    for (auto &ghost:  m_ghosts) { ghost->recv_data(); }

    for (auto &border: m_border) { border->wait_data(); }
    for (auto &ghost:  m_ghosts) { ghost->wait_data(); }
}

void Mesh::sync_adaptation() {
    sync_cells_count();

    bool sync_complete = false;
    while(!sync_complete) {
        sync_flags();

        size_t flags_changed = ghost_flags_changed() ? 1 : 0;
        sync_complete = mpi::max(flags_changed) < 1;

        if (flags_changed) {
            balance_flags();
        }
    }
}

void Mesh::sync_flags() {
    for (auto &border: m_border) { border->send_adaptation_flags(); }
    for (auto &ghost:  m_ghosts) { ghost->recv_adaptation_flags(); }

    for (auto &border: m_border) { border->wait_adaptation_flags(); }
    for (auto &ghost:  m_ghosts) { ghost->wait_adaptation_flags(); }
}

bool Mesh::ghost_flags_changed() {
    for (auto &ghost: m_ghosts) {
        if (ghost->is_adaptation_flags_changed()) {
            return true;
        }
    }
    return false;
}

void Mesh::global_indexing(size_t& idx_begin, size_t& idx_end, size_t& n_size) {
    size_t n_cells = 0;

    for(auto cell: *m_cells) {
        cell->set_index(n_cells);
        ++n_cells;
    }

    vector<size_t> all_sizes = mpi::all_gather(n_cells);
    idx_begin = 0;
    for(int r = 0; r < mpi::rank(); ++r) {
        idx_begin += all_sizes[r];
    }
    idx_end = idx_begin + all_sizes[mpi::rank()];
    n_size = idx_end;
    for(int r = mpi::rank() + 1; r < mpi::size(); ++r) {
        n_size += all_sizes[r];
    }

    for(auto cell: *m_cells) {
        cell->set_index(cell->get_index() + idx_begin);
    }

    for (auto &border: m_border) { border->send_cells_indices(); }
    for (auto &ghost:  m_ghosts) { ghost->recv_cells_indices();  }

    for (auto &border: m_border) { border->wait_cells_indices(); }
    for (auto &ghost:  m_ghosts) { ghost->wait_cells_indices();  }
}

array<size_t, 3> indexing_offsets(array<size_t, 3> sizes) {
    vector<size_t> all_sizes = mpi::all_gather(sizes);

    array<size_t, 3> offsets = { 0, 0, 0 };

    for(int rank = 0; rank < mpi::rank(); ++rank) {
        offsets[0] += all_sizes[rank * 3 + 0];
        offsets[1] += all_sizes[rank * 3 + 1];
        offsets[2] += all_sizes[rank * 3 + 2];
    }
    return offsets;
}

void Mesh::global_indexing(VertexVector& vertices, FaceVector& faces, CellVector& cells) {
    vertices.clear();
    faces.clear();
    cells.clear();

    IndexedObject::index_type n_cells = 0;
    IndexedObject::index_type n_faces = 0;
    IndexedObject::index_type n_verts = 0;

    // Соберем все ячейки
    for(auto cell: *m_cells) {
        cell->set_index(n_cells);
        cells.push_back(cell);
        ++n_cells;
    }

    for(auto cell: *m_fake) {
        cell->set_index(n_cells);
        cells.push_back(cell);
        ++n_cells;
    }

    // Сбросим индексы граней и вершин
    for(const auto& cell: cells) {
        for(auto& face: cell->faces_list()) {
            face->reset_index();
            for(auto& v: face->vertices()) {
                v->reset_index();
            }
        }
    }

    // Установим индексы
    for(auto& cell: cells) {
        for(auto& face: cell->faces_list()) {
            if (face->is_index_undefined()) {
                face->set_index(n_faces);
                faces.push_back(face);
                ++n_faces;
            }

            for(auto& v: face->vertices()) {
                if (v->is_index_undefined()) {
                    v->set_index(n_verts);
                    vertices.push_back(v);
                    ++n_verts;
                }
            }
        }
    }

    // Смещение индексов в глобальной индексации
    auto offsets = indexing_offsets({n_verts, n_faces, n_cells});

    for(auto& vertex: vertices) {
        vertex->set_index(vertex->get_index() + offsets[0]);
    }

    for(auto& face: faces) {
        face->set_index(face->get_index() + offsets[1]);
    }

    for(auto& cell: cells) {
        cell->set_index(cell->get_index() + offsets[2]);
    }

    // синхронизуем индексы
    sync_indices();
}

void Mesh::global_indexing(std::map<IndexedObject::index_type, Vertex::Ptr> &vertices,
                           std::map<IndexedObject::index_type, Face::Ptr> &faces,
                           std::map<IndexedObject::index_type, Cell::Ptr> &cells) {
    VertexVector verts_vector;
    FaceVector faces_vector;
    CellVector cells_vector;

    global_indexing(verts_vector, faces_vector, cells_vector);

    vertices.clear();
    for(const auto& v: verts_vector)
        vertices[v->get_index()] = v;

    faces.clear();
    for(const auto& f: faces_vector)
        faces[f->get_index()] = f;

    cells.clear();
    for(const auto& c: cells_vector)
        cells[c->get_index()] = c;
}

void Mesh::sync_indices() {
    for(auto& ghost:  m_ghosts) { ghost->recv_primitives_indices_count();  }
    for(auto& border: m_border) { border->send_primitives_indices_count(); }

    for(auto& ghost:  m_ghosts) { ghost->wait_primitives_indices_count();  }
    for(auto& border: m_border) { border->wait_primitives_indices_count(); }

    for(auto& ghost:  m_ghosts) { ghost->recv_primitives_indices();  }
    for(auto& border: m_border) { border->send_primitives_indices(); }

    for(auto& ghost:  m_ghosts) { ghost->wait_primitives_indices();  }
    for(auto& border: m_border) { border->wait_primitives_indices(); }
}
