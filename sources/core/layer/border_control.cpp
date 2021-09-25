#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/side.h>
#include <core/cell/cell.h>
#include <core/layer/border_control.h>
#include <core/cell/data_holder.h>
#include <utils/memory/node_list.h>

BorderControl::BorderControl(NodeList::Ptr nodes)
    : LayerControl(std::move(nodes)) {
    m_layer->set_role(ListRole::BORDER);
}

BorderControl::UPtr BorderControl::create(shared_ptr<NodeList> nodes) {
    return makeUnique<BorderControl>(move(nodes));
}

void BorderControl::send_cells_count() {
#ifdef ENABLE_MPI
    m_cells_count = m_layer->size();

    MPI_Isend(&m_cells_count,
              1,
              MPI_SIZE_T,
              m_layer->rank(),
              COUNT_TAG,
              MPI_COMM_WORLD,
              &m_cells_count_request);
#endif
}

void BorderControl::wait_cells_count() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_cells_count_request, MPI_STATUS_IGNORE);
#endif
}

void BorderControl::send_data() {
#ifdef ENABLE_MPI
    auto cell_size = m_layer->begin()->data_holder()->size() + 1;

    m_data.resize(m_layer->size() * cell_size);

    for(size_t i = 0; i < m_sorted_cells.size(); ++i) {
        m_data[cell_size * i] = m_sorted_cells[i]->workload();
        m_sorted_cells[i]->data_holder()->serialize(&m_data[cell_size * i + 1]);

    }

    MPI_Isend(m_data.data(),
              static_cast<int>(m_layer->size() * cell_size),
              MPI_DOUBLE,
              m_layer->rank(),
              DATA_TAG,
              MPI_COMM_WORLD,
              &m_data_request);
#endif
}

void BorderControl::wait_data() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_data_request, MPI_STATUS_IGNORE);
#endif
}

void BorderControl::send_adaptation_flags() {
#ifdef ENABLE_MPI
    m_cells_count = m_layer->size();

    m_adaptation_flags.resize(m_sorted_cells.size());

    for(size_t i = 0; i < m_sorted_cells.size(); ++i) {
        m_adaptation_flags[i] = static_cast<int>(m_sorted_cells[i]->adaptation_flag());
    }

    MPI_Isend(m_adaptation_flags.data(),
              static_cast<int>(m_layer->size()),
              MPI_INT,
              m_layer->rank(),
              ADAPTATION_FLAGS_TAG,
              MPI_COMM_WORLD,
              &m_adaptation_flags_request);
#endif
}

void BorderControl::wait_adaptation_flags() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_adaptation_flags_request, MPI_STATUS_IGNORE);
#endif
}

void BorderControl::send_cells_indices() {
#ifdef ENABLE_MPI
    m_cells_count = m_layer->size();

    m_cells_indices.resize(m_sorted_cells.size());

    for(size_t i = 0; i < m_sorted_cells.size(); ++i) {
        m_cells_indices[i] = static_cast<int>(m_sorted_cells[i]->get_index());
    }

    MPI_Isend(m_cells_indices.data(),
              static_cast<int>(m_layer->size()),
              MPI_SIZE_T,
              m_layer->rank(),
              CELLS_INDICES_TAG,
              MPI_COMM_WORLD,
              &m_cells_indices_request);
#endif
}

void BorderControl::wait_cells_indices() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_cells_indices_request, MPI_STATUS_IGNORE);
#endif
}

void BorderControl::send_primitives_indices_count() {
#ifdef ENABLE_MPI
    m_primitives_indices.clear();

    // На каждую ячейку отправляем
    //        [ID ЯЧЕЙКИ]
    //              [КОЛ-ВО ГРАНЕЙ ПО СТОРОНЕ]
    //              [ID ГРАНЕЙ ПО СТОРОНЕ]
    //                  [ID ВЕРШИН ГРАНЕЙ]
    for(auto& cell: m_sorted_cells) {
        m_primitives_indices.push_back(cell->get_index());


        for (int s = 0; s < 4; ++s) {
            auto faces_list = cell->faces(to_side(s));
            m_primitives_indices.push_back(static_cast<IndexedObject::index_type>(faces_list->size()));

            for (size_t fid = 0; fid < faces_list->size(); ++fid) {
                auto face = faces_list->at(fid);

                m_primitives_indices.push_back(face->get_index());
                m_primitives_indices.push_back(face->vertex(0)->get_index());
                m_primitives_indices.push_back(face->vertex(1)->get_index());
            }
        }
    }

    m_primitives_indices_count = m_primitives_indices.size();

    MPI_Isend(&m_primitives_indices_count,
              1,
              MPI_SIZE_T,
              m_layer->rank(),
              PRIMITIVES_COUNT_TAG,
              MPI_COMM_WORLD,
              &m_primitives_indices_count_request);
#endif
}

void BorderControl::wait_primitives_indices_count() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_primitives_indices_count_request, MPI_STATUS_IGNORE);
#endif
}

void BorderControl::send_primitives_indices() {
#ifdef ENABLE_MPI
    // Данные подготовлены при подсчёте примитивов
    MPI_Isend(m_primitives_indices.data(),
              static_cast<int>(m_primitives_indices.size()),
              MPI_PRIMITIVE_INDEX,
              m_layer->rank(),
              PRIMITIVES_INDICES_TAG,
              MPI_COMM_WORLD,
              &m_primitives_indices_request);
#endif
}

void BorderControl::wait_primitives_indices() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_primitives_indices_request, MPI_STATUS_IGNORE);
#endif
}

void BorderControl::send_double(double value) {
#ifdef ENABLE_MPI
    m_double = value;
    MPI_Isend(&m_double,
              1,
              MPI_DOUBLE,
              m_layer->rank(),
              DOUBLE_TAG,
              MPI_COMM_WORLD,
              &m_double_request);
#endif
}

void BorderControl::wait_double() {
#ifdef ENABLE_MPI
    MPI_Wait(&m_double_request, MPI_STATUS_IGNORE);
#endif
}