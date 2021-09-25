#include <allstd.h>
#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/mesh/mesh.h>
#include <core/layer/ghost_control.h>
#include <core/cell/data_holder.h>
#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/side.h>
#include <core/cell/cell.h>

GhostControl::GhostControl(NodeList::Ptr nodes)
    : LayerControl(std::move(nodes)) {
    m_layer->set_role(ListRole::GHOST);
}

GhostControl::UPtr GhostControl::create(NodeList::Ptr nodes) {
    return makeUnique<GhostControl>(move(nodes));
}

void GhostControl::recv_cells_count() {
#ifdef ENABLE_MPI
    MPI_Irecv(&m_cells_count,
              1,
              MPI_SIZE_T,
              m_layer->rank(),
              COUNT_TAG,
              MPI_COMM_WORLD,
              &m_cells_count_request);
#endif
}

void GhostControl::wait_cells_count() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_cells_count_request, MPI_STATUS_IGNORE);

    if(m_cells_count != m_layer->size()) {
        throw runtime_error("Error! Incorrect ghost size on " + to_string(mpi::rank()) + " get -> " + to_string(m_cells_count) + " expect -> " + to_string(m_layer->size()) + ")!");
    }
#endif
}

void GhostControl::recv_data() {
#ifdef ENABLE_MPI
    auto cell_size = m_layer->begin()->data_holder()->size() + 1;

    m_data.resize(m_layer->size() * cell_size);
    MPI_Irecv(m_data.data(),
              static_cast<int>(m_layer->size() * cell_size),
              MPI_DOUBLE,
              m_layer->rank(),
              DATA_TAG,
              MPI_COMM_WORLD,
              &m_data_request);
#endif
}

void GhostControl::wait_data() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_data_request, MPI_STATUS_IGNORE);

    auto cell_size = m_layer->begin()->data_holder()->size() + 1;
    for(size_t i = 0; i < m_sorted_cells.size(); ++i) {
        m_sorted_cells[i]->reset_workload();
        m_sorted_cells[i]->add_workload(m_data[cell_size * i]);
        m_sorted_cells[i]->data_holder()->deserialize(&m_data[cell_size * i + 1]);
    }
#endif
}

void GhostControl::recv_adaptation_flags() {
#ifdef ENABLE_MPI
    // Если ячеек стало больше - изменим размер буфера флагов
    if (m_layer->size() > m_adaptation_flags.size()) {
        m_adaptation_flags.resize(m_layer->size());
        //m_previous_adaptation_flags.resize(m_layer->size());
    }

    MPI_Irecv(m_adaptation_flags.data(),
              static_cast<int>(m_layer->size()),
              MPI_INT,
              m_layer->rank(),
              ADAPTATION_FLAGS_TAG,
              MPI_COMM_WORLD,
              &m_adaptation_flags_request);
#endif
}

void GhostControl::wait_adaptation_flags() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_adaptation_flags_request, MPI_STATUS_IGNORE);

    m_adaptation_flags_changed = false;

    // Приcвоим флаги и проверим, есть ли в них отличия с последнего приёма
    for(size_t i = 0; i < m_sorted_cells.size(); ++i) {

        /*
        if(m_adaptation_flags[i] != m_previous_adaptation_flags[i]) {
            m_adaptation_flags_changed = true;
        }

        m_previous_adaptation_flags[i] = m_adaptation_flags[i];
         */

        auto recv_flag = static_cast<AdaptationFlag>(m_adaptation_flags[i]);
        if (recv_flag != m_sorted_cells[i]->adaptation_flag()) {
            m_adaptation_flags_changed = true;
        }

        m_sorted_cells[i]->set_adaptation_flag(recv_flag);
    }
#endif
}

bool GhostControl::is_adaptation_flags_changed() {
    return m_adaptation_flags_changed;
}

void GhostControl::recv_cells_indices() {
#ifdef ENABLE_MPI
    // Если ячеек стало больше - изменим размер буфера флагов
    if (m_layer->size() > m_cells_indices.size()) {
        m_cells_indices.resize(m_layer->size());
    }

    MPI_Irecv(m_cells_indices.data(),
              static_cast<int>(m_layer->size()),
              MPI_SIZE_T,
              m_layer->rank(),
              CELLS_INDICES_TAG,
              MPI_COMM_WORLD,
              &m_cells_indices_request);
#endif
}

void GhostControl::wait_cells_indices() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_cells_indices_request, MPI_STATUS_IGNORE);
    for(size_t i = 0; i < m_sorted_cells.size(); ++i) {
        m_sorted_cells[i]->set_index(m_cells_indices[i]);
    }
#endif
}

void GhostControl::recv_primitives_indices_count() {
#ifdef ENABLE_MPI
    MPI_Irecv(&m_primitives_indices_count,
              1,
              MPI_SIZE_T,
              m_layer->rank(),
              PRIMITIVES_COUNT_TAG,
              MPI_COMM_WORLD,
              &m_primitives_indices_count_request);
#endif
}

void GhostControl::wait_primitives_indices_count() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_primitives_indices_count_request, MPI_STATUS_IGNORE);
#endif
}

void GhostControl::recv_primitives_indices() {
#ifdef ENABLE_MPI
    m_primitives_indices.resize(m_primitives_indices_count);

    MPI_Irecv(m_primitives_indices.data(),
              static_cast<int>(m_primitives_indices_count),
              MPI_PRIMITIVE_INDEX,
              m_layer->rank(),
              PRIMITIVES_INDICES_TAG,
              MPI_COMM_WORLD,
              &m_primitives_indices_request);
#endif
}

void GhostControl::wait_primitives_indices() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_primitives_indices_request, MPI_STATUS_IGNORE);

    size_t offset = 0;
    for(auto& cell: m_sorted_cells) {
        cell->set_index(m_primitives_indices[offset++]);

        for(int s = 0; s < 4; ++s) {
            auto side = to_side(s);

            auto faces_list  = cell->faces(side);
            auto faces_count = m_primitives_indices[offset++];

            if (faces_list->size() < faces_count) {
                // Гостовые ячейки не содержат лишних "висящих" вершин,
                // Поэтому может прийти FacesList с бОльшим числом граней
                // Случай реализуется только когда гост ячейка имеет одинарную
                // грань, а получена двойная грань
                auto face = faces_list->single();

                auto face_id1 = m_primitives_indices[offset++];
                auto vert_id1 = m_primitives_indices[offset++];
                offset++; // пропустили среднюю вершину

                auto face_id2 = m_primitives_indices[offset++];
                offset++; // пропустили среднюю вершину
                auto vert_id2 = m_primitives_indices[offset++];

                auto face_id = min(face_id1, face_id2);
                if (face_id < face->get_index()) {
                    face->set_index(face_id);
                }
                if (vert_id1 < face->vertex(0)->get_index()) {
                    face->vertex(0)->set_index(vert_id1);
                }
                if (vert_id2 < face->vertex(1)->get_index()) {
                    face->vertex(1)->set_index(vert_id2);
                }
            }
            else {
                for(size_t fid = 0; fid < faces_count; ++fid) {
                    auto face = faces_list->at(fid);

                    auto face_id = m_primitives_indices[offset++];
                    auto vert_id1 = m_primitives_indices[offset++];
                    auto vert_id2 = m_primitives_indices[offset++];

                    if (face_id < face->get_index()) {
                        face->set_index(face_id);
                    }
                    if (vert_id1 < face->vertex(0)->get_index()) {
                        face->vertex(0)->set_index(vert_id1);
                    }
                    if (vert_id2 < face->vertex(1)->get_index()) {
                        face->vertex(1)->set_index(vert_id2);
                    }
                }
            }
        }
    }
#endif
}

void GhostControl::recv_double() {
#ifdef ENABLE_MPI
    MPI_Irecv(&m_double,
              1,
              MPI_DOUBLE,
              m_layer->rank(),
              DOUBLE_TAG,
              MPI_COMM_WORLD,
              &m_double_request);
#endif
}

double GhostControl::wait_double() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_double_request, MPI_STATUS_IGNORE);
#endif
    return m_double;
}