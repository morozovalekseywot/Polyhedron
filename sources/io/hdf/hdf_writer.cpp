//
// Created by 159-mrv on 10/26/17.
//

#include <hdf5.h>

#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/vertex/vertex.h>
#include <core/cell/cell.h>
#include <core/mesh/mesh.h>
#include <io/hdf/hdf_writer.h>
#include <io/hdf/hdf_utils.h>
#include <utils/memory/node_list.h>
#include "problems/problem.h"

HdfWriter::HdfWriter(const Configuration &config, double init_time, double t_scale)
    : MeshWriter(config, init_time, t_scale) {
}

void HdfWriter::write_head(Problem* problem) {
    HDF::File::create(m_fullname + pt() + to_string(mpi::rank()));

    if (mpi::is_master()) {
        std::ofstream xdmf;
        xdmf.open(m_fullname + ".xmf");

        xdmf << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
        xdmf << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"3.3\">\n";
        xdmf << "  <Domain>\n";
        xdmf << "    <Grid CollectionType=\"Temporal\" GridType=\"Collection\" Name=\"Collection\">\n";

        xdmf.close();
    }
}

void HdfWriter::write_tail() {
    if (mpi::is_master()) {
        std::ofstream xdmf;
        xdmf.open(m_fullname + ".xmf", std::ios_base::app);

        xdmf << "    </Grid>\n";
        xdmf << "  </Domain>\n";
        xdmf << "</Xdmf>\n";

        xdmf.close();
    }
}

#define USE_THREADS 1

void HdfWriter::write(Problem* problem, Mesh *mesh, double time) {
    if (not (time_to_write(time) or (problem && problem->time_to_write()))) {
        return;
    }

    // Сбросить индексы ячеек и вершин
    auto reset_indecies =
            [](const NodeList::Part &cells) {
                for (auto cell: cells) {
                    cell->reset_index();
                    for (auto vertex: cell->ordered_vertices()) {
                        vertex->reset_index();
                    }
                }
            };

#if USE_THREADS
    mesh->map(reset_indecies);
#else
    reset_indecies(mesh->all_chunks());
#endif

    // Занумеруем вершины и посчитаем ячейки
    size_t cells_count = 0;
    size_t points_count = 0;
    for (auto cell: mesh->all_chunks()) {
        cell->set_index(cells_count);
        ++cells_count;
        for (auto &vertex: cell->ordered_vertices()) {
            if (vertex->is_index_undefined()) {
                vertex->set_index(points_count);
                ++points_count;
            }
        }
    }

    // Внесем ячейки и вершины в списки
    vector<float> vertices_table;
    vertices_table.reserve(2 * points_count);
    vector<size_t> cells_table;

    for (auto cell: mesh->all_chunks()) {
        auto ord_vs = cell->ordered_vertices();

        cells_table.push_back(3);
        cells_table.push_back(ord_vs.size());

        for (auto &vertex: ord_vs) {
            size_t v_ind = vertex->get_index();
            cells_table.push_back(v_ind);

            vertices_table[2 * v_ind] = static_cast<float>(vertex->x());
            vertices_table[2 * v_ind + 1] = static_cast<float>(vertex->y());
        }
    }
    size_t cells_table_size = cells_table.size();


    // Открываем файл на запись

    auto hdf5 = HDF::File(m_fullname + pt());

    // Создаем временную группу
    auto time_group = hdf5.add_group("step_" + to_string(m_step));
    time_group.add_attribute("time", H5T_NATIVE_DOUBLE, &time);

    // Создаем таблицы вершнин, ячеек и данных ячеек        
    time_group.add_table("vertices", H5T_NATIVE_FLOAT, {points_count, 2});
    time_group.add_table("cells", H5T_NATIVE_UINT64, {cells_table.size()});

    // Запись геометрии
    time_group.write("vertices", H5T_NATIVE_FLOAT, vertices_table.data());
    time_group.write("cells", H5T_NATIVE_UINT64, cells_table.data());


    // Структура для доступа к данным ячейки
    struct CellData {
        CellData(string p_name, hid_t p_type,
                 std::function<double(const Cell::Ptr &)> p_func) :
                name(move(p_name)), type(p_type), func(std::move(p_func)) {
        }

        string name;
        hid_t type;
        std::function<double(const Cell::Ptr &)> func;
    };

    vector<CellData> cell_data_vector;

    // Добавим параметры солвера
    for(const auto& name: m_variables) {
        cell_data_vector.emplace_back(
                CellData(name, H5T_NATIVE_FLOAT,
                         [problem, &name](const Cell::Ptr &cell) -> double {
                             return problem->get_cell_param(cell, name);
                         })
        );
    }

    // workload
#if 0
    /*
    mesh->map([mesh](const NodeList::Part& cells) {
        for(auto cell:cells) {
            cell->reset_workload();
            cell->add_workload(mesh->cells()->size());
        }
    });
     */

    cell_data_vector.emplace_back(
            CellData("workload", H5T_NATIVE_FLOAT,
            [](const Cell::Ptr&cell) -> double {
                return cell->workload();
            })
    );
#endif

    // owner
#if 0
    cell_data_vector.emplace_back(
            CellData("owner", H5T_NATIVE_INT,
                     [](const Cell::Ptr & cell) -> double {
                         return cell->owner();
                     })
    );
#endif

    // CID
#if 0
    cell_data_vector.emplace_back(
            CellData("cid", H5T_NATIVE_INT,
                     [](const Cell::Ptr&cell) -> double {
                         return cell->id();
                     })
    );
#endif

    // Z
#if 1
    cell_data_vector.emplace_back(
            CellData("cid", H5T_NATIVE_INT,
                     [](const Cell::Ptr&cell) -> double {
                         return cell->z();
                     })
    );
#endif

    // Записываем данные ячеек
    vector<float> fdata(cells_count, 0.0);
    vector<int> idata(cells_count, 0);

    for (auto &cd: cell_data_vector) {
        time_group.add_table(cd.name, cd.type, {cells_count});

#if USE_THREADS
        // Параллельное заполнение массивов данных

        if (cd.type == H5T_NATIVE_FLOAT) {
            mesh->map(
                    [&cd, &fdata](const NodeList::Part &cells) {
                        for (auto cell: cells) {
                            fdata[cell->get_index()] = static_cast<float>(cd.func(cell));
                        }
                    });
            time_group.write(cd.name, cd.type, fdata.data());
        } else {
            mesh->map(
                    [&cd, &idata](const NodeList::Part &cells) {
                        for (auto cell: cells) {
                            idata[cell->get_index()] = static_cast<int>(cd.func(cell));
                        }
                    });
            time_group.write(cd.name, cd.type, idata.data());
        }

#else
        // Последовательное заполнение массивов данных

        if (cd.type == H5T_NATIVE_FLOAT) {
            for (auto cell: mesh->all_chunks()) {
                fdata[cell->get_index()] = static_cast<float>(cd.func(cell));
            }
            time_group.write(cd.name, cd.type, fdata.data());
        } else {
            for (auto cell: mesh->all_chunks()) {
                fdata[cell->get_index()] = static_cast<float>(cd.func(cell));
            }
            time_group.write(cd.name, cd.type, idata.data());
        }
#endif
    }

    time_group.close();
    hdf5.close();

    // Обмениваемся данными о размерах таблиц
    auto ctable_sizes = mpi::all_gather(cells_table_size);
    auto vtable_sizes = mpi::all_gather(points_count);
    auto cells_counts = mpi::all_gather(cells_count);

    if (mpi::is_master()) {
        std::ofstream xdmf;
        xdmf.open(m_fullname + ".xmf", std::ios_base::app);

        xdmf << "<Grid CollectionType=\"None\" GridType=\"Collection\" Name=\"Collection\">";
        xdmf << "<Time Value=\"" << time << "\"/>";

        for (int r = 0; r < mpi::size(); ++r) {
            xdmf << " <Grid Name=\"Grid.Part" << r << "\">";
            xdmf << "<Geometry Origin=\"\" Type=\"XY\">";
            xdmf << "<DataItem DataType=\"Float\" Dimensions=\"" << vtable_sizes[r] << " 2\" "
                 << "Format=\"HDF\" Precision=\"4\">";
            xdmf << m_filename << ".pt" << r << ".hdf5:/step_" << m_step << "/vertices";
            xdmf << "</DataItem>";
            xdmf << "</Geometry>";

            xdmf << "<Topology Dimensions=\"" << cells_counts[r] << "\" Type=\"Mixed\">";
            xdmf << "<DataItem DataType=\"Int\" Dimensions=\"" << ctable_sizes[r] << "\" "
                 << "Format=\"HDF\" Precision=\"8\">";
            xdmf << m_filename << ".pt" << r << ".hdf5:/step_" << m_step << "/cells";
            xdmf << "</DataItem>";
            xdmf << "</Topology>";

            for (const CellData &cd: cell_data_vector) {
                xdmf << "<Attribute Center=\"Cell\" Name=\"" << cd.name << "\" Type=\"Scalar\">";
                xdmf << "<DataItem DataType=\"" << ((cd.type == H5T_NATIVE_FLOAT) ? "Float" : "Int")
                     << "\" Dimensions=\"" << cells_counts[r] << "\" Format=\"HDF\" Precision=\"4\">";
                xdmf << m_filename << ".pt" << r << ".hdf5:/step_" << m_step << "/" << cd.name;
                xdmf << "</DataItem>";
                xdmf << "</Attribute>";
            }
            xdmf << "</Grid>";
        }
        xdmf << "</Grid>\n";

        xdmf.close();
    }

    ++m_step;
}