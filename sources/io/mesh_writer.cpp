#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <io/mesh_writer.h>
#include <io/csv/csv_writer.h>
#include <io/dat/dat_writer.h>
#include <io/vtk/vtk_writer.h>
#include <io/hdf/hdf_writer.h>
#include <problems/problem.h>
#include <utils/math.h>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;

MeshWriter::MeshWriter(const Configuration &config, double init_time, double t_scale) {
    fs::path directory = fs::current_path();
    if (config.exist("output_directory")) {
        fs::path dir = config("output_directory").to_string();
        if (dir.is_relative()) {
            directory /= dir;
        }
        else {
            directory = dir;
        }
    }

    if (!fs::exists(directory) || !fs::is_directory(directory)) {
        fs::create_directories(directory);
    }

    m_filename = config.exist("filename_prefix")  ?
                 config("filename_prefix").to_string() : "mesh";

    m_fullname = (directory/m_filename).string();
    m_tempname = (directory/("." + m_filename)).string();

    if (config.exist("write_frequencies")) {
        for (auto &frequency: config("write_frequencies").to_array()) {
            m_timecodes.push_back(frequency("time").to_double()/t_scale);
            m_frequencies.push_back(frequency("frequency").to_double()/t_scale);
        }
    }
    if (m_timecodes.empty()) {
        m_timecodes   = { std::numeric_limits<double>::max() };
        m_frequencies = { std::numeric_limits<double>::max() };
    }

    double min_freq = *std::min_element(m_frequencies.begin(), m_frequencies.end());

    m_step = 0;
    m_next_write = min_freq * static_cast<unsigned int>(init_time / min_freq);
    m_time_delta = std::numeric_limits<double>::max();

    if (config.exist("variables")) {
        for(auto& param: config("variables").to_array()) {
            m_variables.emplace_back(param.to_string());
        }
    }
}

unique_ptr<MeshWriter> MeshWriter::create(const Configuration &config, double init_time, double t_scale) {
    if (config.exist("output_format")) {
        std::string format = config("output_format").to_string();

        if (format == "CSV") {
            return makeUnique<CsvWriter>(config, init_time, t_scale);
        }
        if (format == "VTK") {
            return makeUnique<VtkWriter>(config, init_time, t_scale);
        }
        if (format == "DAT") {
            return makeUnique<DatWriter>(config, init_time, t_scale);
        }

#ifdef USE_HDF_WRITER
        if (format == "HDF5") {
            return makeUnique<HdfWriter>(config, init_time, t_scale);
        }
#endif

        throw std::runtime_error("Unknown output format '" + format + "'");
    }

    return makeUnique<VtkWriter>(config, init_time, t_scale);
}

bool MeshWriter::time_to_write(double time) {
    bool res = time >= m_next_write;

    if (res) {
        for (size_t i = 0; i < m_timecodes.size(); ++i) {
            if (time >= m_timecodes[i]) {
                double next_timecode = i < m_timecodes.size() - 1 ? m_timecodes[i + 1] : math::inf();
                m_time_delta = min(m_frequencies[i], next_timecode - m_next_write);
            }
        }
        m_next_write += m_time_delta;
    }
    return res;
}

std::string MeshWriter::pt(int rank) {
    if (rank < 0) {
        rank = mpi::rank();
    }
    return mpi::size() > 1 ?  (".pt" + to_string(rank)) : "";
}