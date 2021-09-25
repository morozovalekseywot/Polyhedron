//
// Created by 159-mrv on 11/1/18.
//

#include <control/mpi_wrapper.h>
#include <core/cell/cell.h>
#include <core/mesh/partitioner.h>
#include <core/generator/structured_generator.h>
#include <utils/partitioning/base_splitter_1d.h>
#include <utils/partitioning/base_splitter_2d.h>

Partitioner::Partitioner(const BalancingOptions& options, StructuredGenerator* geometry)
    : m_options(options),
      m_geometry(geometry) {

    // Функция нагрузки
    m_workload = [](Cell::Ref cell) -> double {
        return cell->workload();
    };

    // Определяем Splitter и функцию взятия координат
    if (m_options.cells == CellsSet::BASE) {
        switch (m_options.type) {
            case BalancingType::NONE:
                m_splitter = nullptr;
                break;

            case BalancingType::X:
                m_splitter = makeUnique<BaseSplitter1D>(geometry->global_nx());
                m_location = [this](Cell::Ref cell) -> SplitterFPoint {
                    return {(double) m_geometry->u2xy(cell->z_base())[0], 0.0};
                };
                break;

            case BalancingType::Y:
                m_splitter = makeUnique<BaseSplitter1D>(geometry->global_ny());
                m_location = [this](Cell::Ref cell) -> SplitterFPoint {
                    return {(double) m_geometry->u2xy(cell->z_base())[1], 0.0};
                };
                break;

            case BalancingType::XY:
                m_splitter = makeUnique<BaseSplitter2D>(
                        SplitterIPoint({geometry->global_nx(), geometry->global_ny()})
                );
                m_location = [this](Cell::Ref cell) -> SplitterFPoint {
                    return {
                            (double) m_geometry->u2xy(cell->z_base())[0],
                            (double) m_geometry->u2xy(cell->z_base())[1]
                    };
                };
                break;

            default:
                throw runtime_error("Unknown balancing type");
        }
    } else if (m_options.cells == CellsSet::LEAF) {
        throw runtime_error("Balancing by leaf not supported.");
    } else {
        throw runtime_error("WTF, MAN?");
    }
}

void Partitioner::add_workload(Cell::Ref cell) {
    if (!m_options.use_balancing())
        return;

    m_splitter->add_workload(m_location(cell), m_workload(cell));
}

void Partitioner::sync_data() {
    if (!m_options.use_balancing())
        return;

    m_splitter->sync_data();
}

void Partitioner::split() {
    switch (m_options.type) {
        case BalancingType::NONE:
            return;

        case BalancingType::X:
            m_splitter->split(mpi::size() / m_options.proc_per_y);
            break;

        case BalancingType::Y:
            m_splitter->split(mpi::size() / m_options.proc_per_x);
            break;

        case BalancingType::XY:
            m_splitter->split(mpi::size());
            break;
    }
}

void Partitioner::print_info(const string& tab) const {
    if (!m_options.use_balancing())
        return;

    m_splitter->print_info(tab);
}

int Partitioner::rank(Cell::Ref cell) const {
    if (!m_options.use_balancing()) {
        return 0;
    }

    if (m_options.type == BalancingType::X) {
        size_t part = m_geometry->global_ny() / m_options.proc_per_y + 1;
        if (part * m_options.proc_per_y < m_geometry->global_ny()) {
            ++part;
        }
        int pr_x = m_splitter->rank(m_location(cell));
        int pr_y = int(m_geometry->u2xy(cell->z_base())[1] / part);

        return m_options.proc_per_x * pr_y + pr_x;
    } else if (m_options.type == BalancingType::Y) {
        size_t part = m_geometry->global_nx() / m_options.proc_per_x;
        if (part * m_options.proc_per_x < m_geometry->global_nx()) {
            ++part;
        }
        int pr_x = int(m_geometry->u2xy(cell->z_base())[0] / part);
        int pr_y = m_splitter->rank(m_location(cell));

        return m_options.proc_per_x * pr_y + pr_x;
    } else if (m_options.type == BalancingType::XY) {
        return m_splitter->rank(m_location(cell));
    } else {
        throw runtime_error("Incorrect balancing type.");
    }
}