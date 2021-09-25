#ifndef AMR_BORDER_CONTROL_H
#define AMR_BORDER_CONTROL_H

#include <allstd.h>
#include <control/mpi_wrapper.h>

#include <core/indexed_object.h>
#include <core/layer/layer_control.h>


class BorderControl : public LayerControl {
    using NodeList_SPtr = shared_ptr<NodeList>;

public:
    using UPtr = unique_ptr<BorderControl>;

    explicit BorderControl(NodeList_SPtr nodes);

    ~BorderControl() = default;

    static BorderControl::UPtr create(NodeList_SPtr nodes);

    /** Синхронизация размеров слоев Ghost и Border */
    void send_cells_count();

    void wait_cells_count();


    /** Синхронизация данных */
    void send_data();

    void wait_data();


    /** Синхронизация флагов адаптации */
    void send_adaptation_flags();

    void wait_adaptation_flags();


    /** Синхронизация индексов листовых ячеек */
    void send_cells_indices();

    void wait_cells_indices();


    /** Синхронизация размера буффера примитивов (ячейки, грани и вершины) */
    void send_primitives_indices_count();

    void wait_primitives_indices_count();


    /** Синхронизация индексации ячеек, граней и вершин */
    void send_primitives_indices();

    void wait_primitives_indices();


    /** Отправить число */
    void send_double(double value);

    void wait_double();


private:

    /** Реквесты синхронизации */
#ifdef ENABLE_MPI
    MPI_Request m_cells_count_request;
    MPI_Request m_data_request;
    MPI_Request m_cells_indices_request;
    MPI_Request m_adaptation_flags_request;
    MPI_Request m_primitives_indices_count_request;
    MPI_Request m_primitives_indices_request;
    MPI_Request m_double_request;
#endif


    /** Буфферы для синхронизации */
    size_t m_cells_count;
    vector<double> m_data;
    vector<IndexedObject::index_type> m_cells_indices;
    vector<int> m_adaptation_flags;
    size_t m_primitives_indices_count;
    vector<IndexedObject::index_type> m_primitives_indices;
    double m_double;
};

#endif //AMR_BORDER_CONTROL_H
