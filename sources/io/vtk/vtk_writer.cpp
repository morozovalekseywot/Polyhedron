#include <iostream>
#include <utility>

#include <control/configuration.h>
#include <control/mpi_wrapper.h>

#include <core/cell/data_holder.h>
#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/side.h>
#include <core/cell/cell.h>
#include <core/mesh/mesh.h>

#include <io/vtk/vtk_writer.h>
#include <problems/problem.h>
#include <utils/memory/node_list.h>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <future>

namespace fs = boost::filesystem;

std::vector<std::string> VtkWriter::m_default_variables = {};

bool is_big_endian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

VtkWriter::VtkWriter(const Configuration &config, double init_time, double t_scale)
    : MeshWriter(config, init_time, t_scale) {
    m_default_variables = m_variables;

    fs::remove(m_tempname + ".pvd");
    fs::remove(m_filename + ".pvd");
}

const std::vector<std::string>& VtkWriter::default_params() {
    return m_default_variables;
}

void VtkWriter::write(const NodeList &cells, const string &filename, Problem* problem,
                      const std::vector<std::string>& params) {
    if(cells.empty()) {
        std::cerr << "Empty cells get_list in VTK!" << std::endl;
        return;
    }

    write_binary_xml(problem, cells, filename, params);
}

void VtkWriter::write(Problem* problem, Mesh *mesh, double time) {
    if (not time_to_write(time)) {
        return;
    }

    // write pvd file in new thread
    auto as = std::async(std::launch::async,
                         [this, time]() {
                             write_pvd_step(time);
                         });

    NodeList::Ptr cells = mesh->cells();

    write_binary_xml(problem, *cells);

    as.wait();

    ++m_step;
}

#define HEXAHEDRON_ONLY

//#define WRITE_FLAG
//#define WRITE_OWNER
//#define WRITE_UID
//#define WRITE_INDEX
//#define WRITE_NORMALS
//#define WRITE_NEI_UID
#define WRITE_LEVEL
//#define WRITE_WORKLOAD
//#define WRITE_VOLUME

void VtkWriter::write_binary_header(Problem* problem, const NodeList &cells, std::ofstream &os,
                                    const std::vector<std::string>& params) {
    size_t n_cells = cells.size();
    if (n_cells < 1) {
        throw runtime_error("Error! Trying to write empty cells get_list!");
    }

    bool THREE_DIM = cells.begin()->is_3D();

    // Количество вершин
    size_t N_points = 0;  // Число вершин с учетом кратности
    size_t n_points = 0;  // Число вершин n_points < N_points
    size_t N_faces  = 0;  // Число граней с учетом кратности

    for (auto cell: cells) {
        for (auto& face: cell->faces_list()) {
            for (auto& v: face->vertices()) {
                v->reset_index();
            }
        }
    }

    for (auto cell: cells) {
#ifdef HEXAHEDRON_ONLY
        N_points += (THREE_DIM ? 8 : 4);
#else 
        N_points += cell->ordered_vertices().size();
#endif
        for (auto& face: cell->faces_list()) {
            ++N_faces;
            for(auto& v: face->vertices()) {
                if (v->is_index_undefined()) {
                    v->set_index(n_points++);
                }
            }
        }
    }

    string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";

    os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
    os << "  <UnstructuredGrid>" << '\n';
    os << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells=\"" << n_cells << "\">\n";


    // Points
    size_t offset = 0;
    os << "      <Points>\n";
    os << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
       << offset << "\"/>\n";
    os << "      </Points>\n";
    offset += 3 * n_points * sizeof(double) + sizeof(uint32_t);


    // Cells
    os << "      <Cells>" << '\n';
    os << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += N_points * sizeof(uint64_t) + sizeof(uint32_t);

    os << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\""
       << offset << "\"/>\n";
    offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);

    os << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    offset += n_cells * sizeof(uint8_t) + sizeof(uint32_t);

#ifndef HEXAHEDRON_ONLY
    if (THREE_DIM) {
        os << "        <DataArray type=\"UInt64\" Name=\"faces\" format=\"appended\" offset=\"" << offset << "\"/>\n";
        offset += (n_cells + 5 * N_faces) * sizeof(uint64_t) + sizeof(uint32_t);

        os << "        <DataArray type=\"UInt64\" Name=\"faceoffsets\" format=\"appended\" offset=\"" << offset
           << "\"/>\n";
        offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);
    }
#endif
    os << "      </Cells>\n";


    // CellData
    os << "      <CellData>\n";

    // Данные солвера
    if (problem) {
        for (const auto &name: params) {
            os << "        <DataArray type=\"Float64\" Name=\"" << name
               << "\" format=\"appended\" offset=\"" << offset << "\"/>\n'";
            offset += n_cells * sizeof(double) + sizeof(uint32_t);
        }
    }

#ifdef WRITE_FLAG
    os << "        <DataArray type=\"Int32\" Name=\"flag\" format=\"appended\" offset=\"" << offset << "\"/>" << '\n';
    offset += n_cells * sizeof(int32_t) + sizeof(uint32_t);
#endif

#ifdef WRITE_OWNER
    os << "        <DataArray type=\"Int32\" Name=\"owner\" format=\"appended\" offset=\"" << offset << "\"/>" << '\n';
    offset += n_cells * sizeof(int32_t) + sizeof(uint32_t);
#endif

#ifdef WRITE_UID
    os << "        <DataArray type=\"UInt64\" Name=\"uid\" format=\"appended\" offset=\"" << offset << "\"/>" << '\n';
    offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);
#endif

#ifdef WRITE_INDEX
    os << "        <DataArray type=\"UInt64\" Name=\"index\" format=\"appended\" offset=\"" << offset << "\"/>" << '\n';
    offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);
#endif

#ifdef WRITE_LEVEL
    os << "        <DataArray type=\"UInt32\" Name=\"lvl\" format=\"appended\" offset=\"" << offset << "\"/>" << '\n';
    offset += n_cells * sizeof(uint32_t) + sizeof(uint32_t);
#endif

    static_vector<Side, 6> side_list = { Side::LEFT, Side::RIGHT, Side::BOTTOM, Side::TOP };
    if (THREE_DIM) {
        side_list.push_back(Side::BACK);
        side_list.push_back(Side::FRONT);
    }

#ifdef WRITE_NORMALS
    for(auto side: side_list) {
        os << "        <DataArray type=\"Float64\" Name=\"normal." << to_string(side)
           << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>" << '\n';
        offset += 3 * n_cells * sizeof(double) + sizeof(uint32_t);
    }
#endif

#ifdef WRITE_NEI_UID
    for(auto side: side_list) {
        os << "        <DataArray type=\"UInt64\" Name=\"" << to_string(side)
           << "_nei.uid\" format=\"appended\" offset=\"" << offset << "\"/>" << '\n';
        offset += n_cells * sizeof(uint64_t) + sizeof(uint32_t);
    }
#endif

#ifdef WRITE_WORKLOAD
    os << "        <DataArray type=\"Float64\" Name=\"workload\" format=\"appended\" offset=\"" << offset << "\"/>" << '\n';
    offset += n_cells * sizeof(double) + sizeof(uint32_t);
#endif

#ifdef WRITE_VOLUME
    os << "        <DataArray type=\"Float64\" Name=\"volume\" format=\"appended\" offset=\"" << offset << "\"/>" << '\n';
    offset += n_cells * sizeof(double) + sizeof(uint32_t);
#endif

    os << "      </CellData>\n";
    os << "    </Piece>\n";
    os << "  </UnstructuredGrid>\n";

    // AppendedData
    os << "  <AppendedData encoding=\"raw\">\n";
    os << "_";
}

void VtkWriter::write_binary_data(Problem* problem, const NodeList &cells, std::ofstream &os,
                                  const std::vector<std::string>& params) {
    size_t n_cells = cells.size();

    bool THREE_DIM = cells.begin()->is_3D();

    // Количество вершин
    size_t N_points = 0;  // Число вершин с учетом кратности
    size_t n_points = 0;  // Число вершин n_points < N_points
    size_t N_faces = 0;  // Число граней с учетом кратности

    for (auto cell: cells) {
        for (auto &face: cell->faces_list()) {
            for (auto &v: face->vertices()) {
                v->reset_index();
            }
        }
    }

    for (auto cell: cells) {
#ifdef HEXAHEDRON_ONLY
        N_points += (THREE_DIM ? 8 : 4);
#else
        N_points += cell->ordered_vertices().size();
#endif
        for (auto &face: cell->faces_list()) {
            ++N_faces;
            for (auto &v: face->vertices()) {
                if (v->is_index_undefined()) {
                    v->set_index(n_points++);
                }
            }
        }
    }

    // PointsCoords
    uint32_t data_size = static_cast<uint32_t>(3 * n_points * sizeof(double));
    os.write((char *) &data_size, sizeof(uint32_t));

    for (auto cell: cells) {
        for (auto &face: cell->faces_list()) {
            for (auto &v: face->vertices()) {
                v->reset_index();
            }
        }
    }
    size_t counter = 0;
    for (auto cell: cells) {
        for (auto &face: cell->faces_list()) {
            for (auto &v: face->vertices()) {
                if (v->is_index_undefined()) {
                    Vector3d coord = v->v();
                    os.write((char *) &coord, 3 * sizeof(double));
                    v->set_index(counter++);
                }
            }
        }
    }

    // Cells
    // Connectivity
    data_size = static_cast<uint32_t>(N_points * sizeof(uint64_t));
    os.write((char *) &data_size, sizeof(uint32_t));

    for (auto cell: cells) {
#ifdef HEXAHEDRON_ONLY
        vector<Vertex::Ptr> vertices;
        if (THREE_DIM) {
            array<Vertex::Ptr, 4> lv;
            array<Vertex::Ptr, 4> rv;

            auto left = cell->faces(Side::LEFT);
            for (int i = 0; i < 4; ++i) {
                uint k = left->is_simple() ? 0 : (uint) i;
                lv[i] = left->at(k)->vertex(i);
            }
            auto right = cell->faces(Side::RIGHT);
            for (int i = 0; i < 4; ++i) {
                uint k = right->is_simple() ? 0 : (uint) i;
                rv[i] = right->at(k)->vertex(i);
            }

            // Нормали граней направлены в разные стороны
            if (cell->is_correct_face(Side::LEFT) ==
                    cell->is_correct_face(Side::RIGHT)) {
                std::swap(rv[1], rv[3]);
            }

            int zero_idx = 0;
            double min_dist2 = math::inf();
            for (int i = 0; i < 4; ++i) {
                double dist2 = (rv[i]->v() - lv[0]->v()).squaredNorm();
                if (dist2 < min_dist2) {
                    min_dist2 = dist2;
                    zero_idx = i;
                }
            }

            vertices.resize(8);
            for(int i = 0; i < 4; ++i) {
                vertices[i] = lv[i];
                vertices[4 + i] = rv[(4 + i + zero_idx) % 4];
            }
        }
        else {
            auto left = cell->faces(Side::LEFT);
            vertices.push_back(left->vertex(0));
            vertices.push_back(left->vertex(1));

            auto right = cell->faces(Side::RIGHT);
            vertices.push_back(right->vertex(0));
            vertices.push_back(right->vertex(1));

            if (cell->is_correct_face(Side::LEFT) != cell->is_correct_face(Side::RIGHT)) {
                std::swap(vertices[2], vertices[3]);
            }
        }
#else
        auto vertices = cell->ordered_vertices();
#endif
        
        vector<uint64_t> tmps(vertices.size());
        for (size_t k = 0; k < vertices.size(); ++k) {
            tmps[k] = vertices[k]->get_index();
        }

        os.write((char *) tmps.data(), sizeof(uint64_t) * vertices.size());
    }

    // Offsets
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
    os.write((char *) &data_size, sizeof(uint32_t));

    uint64_t offset = 0;
    for (auto cell: cells) {
#ifdef HEXAHEDRON_ONLY
        if (THREE_DIM) {
            offset += 8;
        }
        else {
            offset += 4;
        }
#else
        offset += cell->ordered_vertices().size();
#endif
        os.write((char *) &(offset), sizeof(uint64_t));
    }

    // Types
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint8_t));
    os.write((char *) &data_size, sizeof(uint32_t));

    uint8_t type;
#ifdef HEXAHEDRON_ONLY
    type = static_cast<uint8_t>(THREE_DIM ? 12 : 9);
#else
    type = static_cast<uint8_t>(THREE_DIM ? 42 : 7);
#endif
    for (size_t i = 0; i < n_cells; ++i) {
        os.write((char *) &type, sizeof(uint8_t));
    }

#ifndef HEXAHEDRON_ONLY
    if (THREE_DIM) {
        // Faces
        data_size = static_cast<uint32_t>((n_cells + 5 * N_faces) * sizeof(uint64_t));
        os.write((char *) &data_size, sizeof(uint32_t));
        for (auto cell: cells) {
            auto faces = cell->faces();
            uint64_t nf = faces.size();
            os.write((char *) &nf, sizeof(uint64_t));
            for (auto &face: faces) {
                array<uint64_t, 5> fi = {
                        4,
                        face->vertices()[0]->get_index(),
                        face->vertices()[1]->get_index(),
                        face->vertices()[2]->get_index(),
                        face->vertices()[3]->get_index()
                };
                os.write((char *) fi.data(), 5 * sizeof(uint64_t));
            }
        }

        // FaceOffsets
        data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
        os.write((char *) &data_size, sizeof(uint32_t));
        offset = 0;
        for (auto cell: cells) {
            offset += (1 + 5 * cell->faces().size());
            os.write((char *) &(offset), sizeof(uint64_t));
        }
    }
#endif


    // CellData
    if (problem) {
        for (const auto &name: params) {
            data_size = static_cast<uint32_t>(n_cells * sizeof(double));
            os.write((char *) &data_size, sizeof(uint32_t));

            for (auto cell: cells) {
                double value = 0.0;
                if (cell->data_holder()->size() < 1) {
                    value = -777.0;
                } else {
                    value = problem->get_cell_param(cell, name);
                }
                os.write((char *) &value, sizeof(double));
            }
        }
    }


#ifdef WRITE_FLAG
    data_size = static_cast<uint32_t>(n_cells * sizeof(int32_t));
    os.write((char *) &data_size, sizeof(uint32_t));
    for (auto cell: cells) {
        int32_t value = static_cast<int32_t>(cell->adaptation_flag());
        os.write((char *) &(value), sizeof(int32_t));
    }
#endif

#ifdef WRITE_OWNER
    data_size = static_cast<uint32_t>(n_cells * sizeof(int32_t));
    os.write((char*)&data_size, sizeof(int32_t));
    for(auto cell: cells) {
        int  value = cell->owner();
        os.write((char*)&(value), sizeof(int32_t));
    }
#endif

#ifdef WRITE_UID
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
    os.write((char*)&data_size, sizeof(uint32_t));
    for(auto cell: cells) {
        uint64_t value = cell->id();
        os.write((char*)&(value), sizeof(uint64_t));
    }
#endif

#ifdef WRITE_INDEX
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));
    os.write((char*)&data_size, sizeof(uint32_t));
    for(auto cell: cells) {
        uint64_t value = cell->get_index();
        os.write((char*)&(value), sizeof(uint64_t));
    }
#endif

#ifdef WRITE_LEVEL
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint32_t));
    os.write((char*)&data_size, sizeof(uint32_t));
    for(auto cell: cells) {
        uint32_t value = static_cast<uint32_t>(cell->level());
        os.write((char*)&(value), sizeof(uint32_t));
    }
#endif

    static_vector<Side, 6> side_list = {Side::LEFT, Side::RIGHT, Side::BOTTOM, Side::TOP};
    if (THREE_DIM) {
        side_list.push_back(Side::BACK);
        side_list.push_back(Side::FRONT);
    }

#ifdef WRITE_NORMALS
    data_size = static_cast<uint32_t>(3 * n_cells * sizeof(double));
    for (auto side: side_list) {
        os.write((char *) &data_size, sizeof(uint32_t));
        for (auto cell: cells) {
            auto face = cell->faces(side)->single();
            Vector3d coord = { NAN, NAN, NAN };
            //Vector3d coord = face->center() - cell->center();
            //coord.normalize();
            if (face->neighbor_exist(cell)) {
                auto nei = face->neighbor(cell);
                coord = face->normal(nei);
            }
            //Vector3d coord = face->vertex(0)->v();
            //coord = (cell->is_correct_face(side) ? 1.0 : -1.0) * face->outward_normal();
            //if (face->neighbor_exist(cell)) {
            //    coord = face->neighbor(cell)->center() - cell->center();
            //    coord.normalize();
            //}
            //coord = face->outward_normal();
            os.write((char *) coord.data(), 3 * sizeof(double));
        }
    }
#endif

#ifdef WRITE_NEI_UID
    data_size = static_cast<uint32_t>(n_cells * sizeof(uint64_t));

    for(auto side: side_list) {
        os.write((char *) &data_size, sizeof(uint32_t));
        for (auto cell: cells) {
            auto face = cell->faces(side)->single();
            uint64_t nei_uid = 0;
            if (face->neighbor_exist(cell)) {
                nei_uid = face->neighbor(cell)->id();
            }
            os.write((char *)&nei_uid, sizeof(uint64_t));
        }
    }
#endif

#ifdef WRITE_WORKLOAD
    data_size = static_cast<uint32_t>(cells_count * sizeof(double));

    os.write((char *) &data_size, sizeof(uint32_t));
    for (auto cell: cells) {
        auto batya = Cell::origin(cell);
        double w = batya->workload();
        os.write((char *) &w, sizeof(double));
    }
#endif

#ifdef WRITE_VOLUME
    data_size = static_cast<uint32_t>(n_cells * sizeof(double));

    os.write((char *) &data_size, sizeof(uint32_t));
    for (auto cell: cells) {
        double w = cell->volume();
        os.write((char *) &w, sizeof(double));
    }
#endif
    // CellData end
}

void VtkWriter::write_binary_footer(std::ofstream &os) {
    os << "  </AppendedData>\n";
    os << "</VTKFile>\n";
}

void VtkWriter::write_binary_xml(Problem* problem, const NodeList &cells, const string &filename,
                                 const std::vector<std::string>& params) {
    std::ofstream os;
    os.open(filename, std::ios::out | std::ios::binary);

    write_binary_header(problem, cells, os, params);
    write_binary_data  (problem, cells, os, params);
    write_binary_footer(os);

    os.close();
}

void VtkWriter::write_binary_xml(Problem* problem, const NodeList &cells) {
    string filename = m_fullname +
            "_" + to_string(m_step) + pt() + ".vtu";

    std::ofstream os;
    os.open(filename, std::ios::out | std::ios::binary);

    write_binary_header(problem, cells, os, m_variables);
    write_binary_data  (problem, cells, os, m_variables);
    write_binary_footer(os);

    os.close();
}

void VtkWriter::write_head(Problem* problem) { }

void VtkWriter::write_tail() { }

void VtkWriter::write_pvd_step(double time) {
    if (!mpi::is_master()) {
        return;
    }
    string fullname = m_fullname + ".pvd";
    string tempname = m_tempname + ".pvd";

    // Create m_tempname.pvd and m_fullname.pvd if they don't exist
    if (!fs::exists(fullname) || !fs::exists(tempname)) {
        std::ofstream os(tempname, std::ios_base::out);
        string byteord = is_big_endian() ? "BigEndian" : "LittleEndian";
        os << "<?xml version=\"1.0\"?>\n";
        os << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" + byteord + "\">\n";
        os << "\t<Collection>\n";
        os.close();

        fs::copy_file(tempname, fullname, fs::copy_option::overwrite_if_exists);
        os.open(fullname, std::ios_base::app);
        os << "\t</Collection>\n";
        os << "</VTKFile>\n";
        os.close();
    }

    std::ofstream os(tempname, std::ios_base::app);
    os << std::scientific << std::setprecision(4);
    for (size_t rank = 0; rank < static_cast<size_t>(mpi::size()); ++rank) {
        os << "\t\t<DataSet timestep=\"" << time <<
           "\" part=\"" + to_string(rank) +
           "\" file=\"" + m_filename + "_" + to_string(m_step) + pt(rank) + ".vtu\"/>\n";
    }
    os.close();

    fs::copy_file(tempname, fullname, fs::copy_option::overwrite_if_exists);
    os.open(fullname, std::ios_base::app);
    os << "\t</Collection>\n";
    os << "</VTKFile>\n";
    os.close();
}
