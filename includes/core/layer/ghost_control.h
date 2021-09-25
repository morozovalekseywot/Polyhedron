#ifndef AMR_GHOST_CONTROL_H
#define AMR_GHOST_CONTROL_H

#include <allstd.h>
#include <control/mpi_wrapper.h>

#include <core/indexed_object.h>
#include <core/layer/layer_control.h>

class NodeList;

class GhostControl : public LayerControl {

public:
    using UPtr = unique_ptr<GhostControl>;
    using URef = const unique_ptr<GhostControl> &;

    explicit GhostControl(shared_ptr<NodeList> nodes);

    ~GhostControl() = default;

    static GhostControl::UPtr create(shared_ptr<NodeList> nodes);


    /** Синхронизация размеров слоев Ghost и Border */
    void recv_cells_count();

    void wait_cells_count();


    /** Синхронизация данных */
    void recv_data();

    void wait_data();


    /** Синхронизация флагов адаптации */
    void recv_adaptation_flags();

    void wait_adaptation_flags();

    bool is_adaptation_flags_changed();


    /** Синхронизация индексов листовых ячеек */
    void recv_cells_indices();

    void wait_cells_indices();


    /** Синхронизация количества примитивов */
    void recv_primitives_indices_count();

    void wait_primitives_indices_count();


    /** Синхронизация индексации граней и вершин */

    void recv_primitives_indices();

    void wait_primitives_indices();


    /** Получить число */
    void recv_double();

    double wait_double();


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
    //vector<int> m_previous_adaptation_flags;
    bool m_adaptation_flags_changed;
    size_t m_primitives_indices_count;
    vector<IndexedObject::index_type> m_primitives_indices;
    double m_double;
};

#endif //AMR_GHOST_CONTROL_H
