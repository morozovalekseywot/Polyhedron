#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/cell/cell.h>
#include <core/face/face.h>
#include <core/vertex/vertex.h>
#include <core/mesh/mesh.h>
#include <io/csv/csv_writer.h>
#include <problems/problem.h>


CsvWriter::CsvWriter(const Configuration &config, double init_time, double t_scale)
        : MeshWriter(config, init_time, t_scale) {
}

void CsvWriter::write_head(Problem* problem) { }

void CsvWriter::write_tail() { }

void CsvWriter::write(Problem* problem, Mesh *mesh, double time) {
    if (not time_to_write(time)) {
        return;
    }

    string filename = m_fullname + "_" + to_string(m_step) + pt() + ".csv";

    std::ofstream os;
    os.open(filename, std::ios::out | std::ios::binary);

    auto cells = mesh->cells();

    os << std::scientific << std::setprecision(6);
    os << "# n_cells, time\n";
    os << cells->size() << "," << time << "\n";
    os << "# x_min,x_max,y_min,y_max,z_min,z_max";
    for (auto& p: m_variables) {
        os << "," << p;
    }
    os << "\n";

    for(auto cell: *cells) {
        double x_min = std::numeric_limits<double>::max();
        double y_min = std::numeric_limits<double>::max();
        double z_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::min();
        double y_max = std::numeric_limits<double>::min();
        double z_max = std::numeric_limits<double>::min();
        for (auto& face: cell->faces_list()) {
            for(auto& v: face->vertices()) {
                x_min = std::min(x_min, v->x());
                y_min = std::min(y_min, v->y());
                z_min = std::min(z_min, v->z());
                x_max = std::max(x_max, v->x());
                y_max = std::max(y_max, v->y());
                z_max = std::max(z_max, v->z());
            }
        }

        os << x_min << "," << x_max << "," << y_min << "," << y_max << "," << z_min << "," << z_max;
        for (auto& p: m_variables) {
            os << "," << problem->get_cell_param(cell, p);
        }

        os << "\n";
    }

    os.close();

    ++m_step;
}