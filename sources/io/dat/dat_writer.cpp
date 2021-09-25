#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/cell/cell.h>
#include <core/face/face.h>
#include <core/vertex/vertex.h>
#include <core/mesh/mesh.h>
#include <io/dat/dat_writer.h>
#include <problems/problem.h>


DatWriter::DatWriter(const Configuration &config, double init_time, double t_scale)
        : MeshWriter(config, init_time, t_scale) {
}

void DatWriter::write_head(Problem* problem) {
    if (!mpi::is_master()) {
        return;
    }

    string filename = m_fullname + ".dat";
    std::ofstream os;
    os.open(filename, std::ios::out | std::ios::trunc);

    os << "# time";
    for (auto& name: m_variables) {
        os << "\t" << name;
    }
    os << "\n";

    os.close();
}

void DatWriter::write_tail() { }

void DatWriter::write(Problem* problem, Mesh *mesh, double time) {
    if (!mpi::is_master()) {
        return;
    }
    if (not time_to_write(time)) {
        return;
    }

    string filename = m_fullname + ".dat";

    std::ofstream os;
    os.open(filename, std::ios::app);

    os << std::scientific << std::setprecision(6) << time;
    for (auto& name: m_variables) {
        os << "\t" << problem->get_integral_param(name);
    }
    os << "\n";
    os.close();

    ++m_step;
}