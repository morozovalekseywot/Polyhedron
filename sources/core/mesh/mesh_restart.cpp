//
// Created by 159-mrv on 4/20/18.
//

#include <control/configuration.h>
#include <control/mpi_wrapper.h>

#include <core/vertex/vertex.h>
#include <core/face/face.h>
#include <core/face/faces_list.h>
#include <core/cell/cell.h>
#include <core/cell/data_holder.h>

#include <core/layer/border_control.h>
#include <core/layer/ghost_control.h>
#include <core/mesh/mesh.h>
#include <core/generator/structured_decomposition.h>
#include <core/generator/sector_generator.h>

#include <io/vtk/vtk_writer.h>
#include <problems/problem.h>

#include <core/mesh/partitioner.h>


static const string vertices_str = "vertices";
static const string faces_str = "faces";
static const string cells_str = "cells";
static const string data_str = "data";

static const std::set<string> sections = {
        vertices_str,
        faces_str,
        cells_str,
        data_str
};

void Mesh::checkpoint(string path, double time, bool binary) {
    // ========================================================================
    // Создадим основной файл чекпоинта и определим имя файла процесса
    // ========================================================================

    if (mpi::is_master()) {
        auto ind = path.find_last_of('/');
        string prefix = (ind != path.size()) ? path.substr(ind + 1) : path;

        std::ofstream os(path + ".chp", std::ios_base::out);

        os << "time: " << std::setprecision(30) << time << "\n";
        os << "format: " << ( binary ? "binary" : "ascii" ) << "\n";
        for (int r = 0; r < mpi::size(); ++r) {
            os << prefix << ".pt" << r << ".chp\n";
        }
    }
    path += ".pt" + mpi::srank() + ".chp";


    // ========================================================================
    // Индексируем примитивы и собираем их в вектора
    // ========================================================================

    vector<Vertex::Ptr> verts;
    vector<Face::Ptr>   faces;
    vector<Cell::Ptr>   cells;

    global_indexing(verts, faces, cells);


    // ========================================================================
    // Сохраняем файл
    // ========================================================================


    if (binary) {
        std::ofstream file(path, std::ios_base::out | std::ios_base::binary);

        file.write(vertices_str.data(), vertices_str.size());
        {
            size_t n_verts = verts.size();
            vector<ClearVertex> clear_verts(n_verts);
            for (size_t i = 0; i < n_verts; ++i) {
                clear_verts[i] = ClearVertex(verts[i]);
            }
            file.write(reinterpret_cast<const char *>(&n_verts), sizeof(n_verts));
            file.write(reinterpret_cast<const char *>(clear_verts.data()),
                       n_verts * sizeof(ClearVertex));
        }

        file.write(faces_str.data(), faces_str.size());
        {
            size_t n_faces = faces.size();
            vector<ClearFace> clear_faces(n_faces);
            for(size_t i = 0; i < n_faces; ++i) {
                clear_faces[i] = ClearFace(faces[i]);
            }
            file.write(reinterpret_cast<const char *>(&n_faces), sizeof(n_faces));
            file.write(reinterpret_cast<const char *>(clear_faces.data()),
                       n_faces * sizeof(ClearFace));
        }

        file.write(cells_str.data(), cells_str.size());
        {
            size_t n_cells = cells.size();
            vector<ClearCell> clear_cells(n_cells);
            for(size_t i = 0; i < n_cells; ++i) {
                clear_cells[i] = ClearCell(cells[i], m_fake->contain(cells[i]));
            }
            file.write(reinterpret_cast<const char *>(&n_cells), sizeof(n_cells));
            file.write(reinterpret_cast<const char *>(clear_cells.data()),
                       n_cells * sizeof(ClearCell));
        }

        file.write(data_str.data(), data_str.size());
        {
            size_t n_data = cells.size();
            file.write(reinterpret_cast<const char *>(&n_data), sizeof(n_data));

            for (auto &cell: cells) {
                size_t data_size = cell->data_holder()->size();

                file.write(reinterpret_cast<const char *>(&data_size), sizeof(data_size));
                file.write(reinterpret_cast<const char *>(cell->data_holder()->buffer()),
                           data_size * sizeof(double));
            }
        }

        file.close();
    }
    else {
        std::ofstream file(path, std::ios_base::out);

        /*
        file << vertices_str << "\n\n";
        for (auto &vertex: verts) {
            file << ClearVertex(vertex) << "\n";
        }

        file << faces_str << "\n\n";
        for (auto &face: faces) {
            file << ClearFace(face) << "\n";
        }

        file << cells_str << "\n\n";
        for (auto &cell: cells) {
            file << ClearCell(cell, m_fake->contain(cell)) << "\n";
        }

        file << data_str << "\n\n";
        for (auto &cell: cells) {
            auto data_size = cell->data_holder()->size();
            double *data = cell->data_holder()->buffer();

            file << data_size;
            for (size_t i = 0; i < data_size; ++i) {
                file << " " << std::setprecision(30) << data[i];
            }
            file << "\n";
        }
         */

        for(auto cell: *m_cells) {
            size_t data_size = cell->data_holder()->size();
            double* data = cell->data();
            file << cell->get_index() << " " << cell->level() << " "
                 << cell->z() << " "
                 << data_size << " " << setprecision(20);
            for(size_t i = 0; i < data_size; ++i) {
                file << data[i] << " ";
            }
            file << "\n";
        }

        file.close();
    }
}

double Mesh::checkpoint_converter(string full_path) {
    // ========================================================================
    // Чтение основного файла чекпоинта
    // ========================================================================

    // Рестарт
    mpi::cout << "  Start calculation from checkpoint '" << full_path << "'.\n";

    double restart_time;
    bool binary;
    vector<string> partition_files;
    {
        auto max_stream_size = std::numeric_limits<std::streamsize>::max();
        string line;

        string directory;
        auto ind = full_path.find_last_of('/');
        if (ind != full_path.size()) {
            directory = full_path.substr(0, ind + 1);
        }

        std::ifstream main_file(full_path, std::ios::in);
        main_file.ignore(max_stream_size, ' ');
        main_file >> restart_time;
        main_file.ignore(max_stream_size, '\n');

        main_file.ignore(max_stream_size, ' ');
        std::getline(main_file, line);
        if (line == "binary") {
            binary = true;
        } else if (line == "ascii") {
            binary = false;
        } else {
            throw runtime_error("Checkpoint syntax error. Wrong checkpoint format '" + line + "', "
                    "avaliable: 'binary'/'ascii'.");
        }

        while (std::getline(main_file, line)) {
            partition_files.push_back(directory + line);
        }
        main_file.close();
    }


    // ========================================================================
    // Считываем примитивы из чекпоинтов в вектора
    // ========================================================================
    vector<ClearVertex> clear_verts;
    vector<ClearFace> clear_faces;
    vector<ClearCell> clear_cells;
    vector<double> clear_data;

    for (auto &filename: partition_files) {
        mpi::cout << "    Load: " << filename << "\n";

        if (binary) {
            std::ifstream file(filename, std::ios_base::binary);

            if (!file) {
                throw runtime_error("Can't open checkpoint file " + filename);
            }

            file.ignore(vertices_str.size());
            {
                size_t n_verts;
                size_t n_verts_old = clear_verts.size();
                file.read(reinterpret_cast<char *>(&n_verts), sizeof(n_verts));

                clear_verts.resize(n_verts_old + n_verts);
                file.read(reinterpret_cast<char *>(&clear_verts[n_verts_old]),
                          n_verts * sizeof(ClearVertex));
            }

            file.ignore(faces_str.size());
            {
                size_t n_faces;
                size_t n_faces_old = clear_faces.size();
                file.read(reinterpret_cast<char *>(&n_faces), sizeof(n_faces));

                clear_faces.resize(n_faces_old + n_faces);
                file.read(reinterpret_cast<char *>(&clear_faces[n_faces_old]),
                          n_faces * sizeof(ClearFace));
            }

            file.ignore(cells_str.size());
            {
                size_t n_cells;
                size_t n_cells_old = clear_cells.size();
                file.read(reinterpret_cast<char *>(&n_cells), sizeof(n_cells));

                clear_cells.resize(n_cells_old + n_cells);
                file.read(reinterpret_cast<char *>(&clear_cells[n_cells_old]),
                          n_cells * sizeof(ClearCell));
            }

            file.ignore(data_str.size());
            {
                size_t n_cells;
                file.read(reinterpret_cast<char *>(&n_cells), sizeof(n_cells));

                for (size_t i = 0; i < n_cells; ++i) {
                    size_t data_size;
                    file.read(reinterpret_cast<char *>(&data_size), sizeof(data_size));
                    vector<double> temp(data_size);
                    file.read(reinterpret_cast<char *>(temp.data()),
                              data_size * sizeof(double));

                    for (auto d: temp) {
                        clear_data.push_back(d);
                    }
                }
            }

            file.close();
        } else {
            std::ifstream file(filename);

            string line;
            string section;
            while (std::getline(file, line)) {
                if (!line.empty()) {
                    if (sections.count(line) > 0) {
                        section = line;
                        continue;
                    }

                    if (section == vertices_str) {
                        clear_verts.emplace_back(ClearVertex(line));

                    } else if (section == faces_str) {
                        clear_faces.emplace_back(ClearFace(line));

                    } else if (section == cells_str) {
                        clear_cells.emplace_back(ClearCell(line));

                    } else if (section == data_str) {
                        std::istringstream ss(line);
                        size_t data_size;
                        ss >> data_size;

                        for (size_t i = 0; i < data_size; ++i) {
                            double data;
                            ss >> data;
                            clear_data.push_back(data);
                        }
                    }
                }
            }

            file.close();
        }

        // Считали новый файл, дай конвертну в новый формат
        ofstream file(filename + ".new", std::ios::out | std::ios::binary);
        size_t data_size = clear_data.size() / clear_cells.size();
        std::cout << "cd.size(): " << clear_data.size() << "; cc.size(): " << clear_cells.size() << "\n";
        for(size_t i = 0; i < clear_cells.size(); ++i) {
            if (clear_cells[i].fake) {
                continue;
            }

            file << clear_cells[i].index << " " << clear_cells[i].level << " "
                 << clear_cells[i].z << " "
                 << data_size << " " << setprecision(20);
            for(size_t j = 0; j < data_size; ++j) {
                file << clear_data[data_size * i + j] << " ";
            }
            file << "\n";
        }
        file.close();
        clear_cells.clear();
        clear_verts.clear();
        clear_faces.clear();
        clear_data.clear();
    }

    mpi::barrier();
    mpi::cout << "    Success load\n";
    throw runtime_error("Enough");
}

/// Рестарт будем делать на основе созданной базовой сетки
/// методом последовательной адаптации
double Mesh::restart_from_checkpoint(string full_path) {
    // ========================================================================
    // Чтение основного файла чекпоинта
    // ========================================================================

    // Рестарт
    mpi::cout << "  Start calculation from checkpoint '" << full_path << "'.\n";

    double restart_time;
    bool binary;
    vector<string> partition_files;
    {
        auto max_stream_size = std::numeric_limits<std::streamsize>::max();
        string line;

        string directory;
        auto ind = full_path.find_last_of('/');
        if (ind != full_path.size()) {
            directory = full_path.substr(0, ind + 1);
        }

        std::ifstream main_file(full_path, std::ios::in);
        main_file.ignore(max_stream_size, ' ');
        main_file >> restart_time;
        main_file.ignore(max_stream_size, '\n');

        main_file.ignore(max_stream_size, ' ');
        std::getline(main_file, line);
        if (line == "binary") {
            binary = true;
        } else if (line == "ascii") {
            binary = false;
        } else {
            throw runtime_error("Checkpoint syntax error. Wrong checkpoint format '" + line + "', "
                    "avaliable: 'binary'/'ascii'.");
        }

        while (std::getline(main_file, line)) {
            partition_files.push_back(directory + line);
        }
        main_file.close();
    }

    size_t data_size;
    vector<double> temp;
    // cells_to_split[l] - ячейки l-го уровня, которые должны разбиться до l+1-го
    vector<std::set<size_t>> cells_to_split(m_max_level);

    for(auto& part_filename: partition_files) {
        std::cout << "PART FN: " << part_filename << "\n";
        ifstream file;
        if (binary) {
            file.open(part_filename, std::ios::in | std::ios::binary);
        }
        else {
            file.open(part_filename, std::ios::in);
        }

        while(!file.eof()) {
            size_t z_top, level;
            IndexedObject::index_type idx;

            if (binary) {
                file.read((char*)&idx, sizeof(idx));
                file.read((char*)&level, sizeof(level));
                file.read((char*)&z_top, sizeof(z_top));
                file.read((char*)&data_size, sizeof(data_size));
                temp.resize(data_size);
                file.read((char*)temp.data(), data_size * sizeof(double));
            }
            else {
                file >> idx >> level >> z_top >> data_size;
                file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }

            size_t z = z_top;
            for (int l = (int)level - 1; l >= 0; --l) {
                z /= 4;
                cells_to_split[l].insert(z);
            }
        }
        file.close();
    }

    extend_cell_data(data_size);
    mpi::cout << "    Refine to\n";
    for(uint l = 0; l < m_max_level; ++l) {
        mpi::cout << "      level " << l + 1 << "\n";

        // set_flags
        for(auto cell: *m_cells) {
            reset_adaptation_flag(cell);
        }
        for(auto cell: *m_cells) {
            if (cell->level() != l) {
                continue;
            }
            if (cells_to_split[l].count(cell->z()) > 0) {
                set_need_to_split(cell);
            }
        }

        sync_adaptation();

        apply_flags(nullptr);

        for(auto& border: m_border) { border->rebuild(); }
        for(auto& ghost:  m_ghosts) { ghost->rebuild();  }

        load_balance();

        m_cells->split(m_n_chunks);
    }
    cells_to_split.clear();
    cells_to_split.shrink_to_fit();

    m_cells->shuffle();

    m_cells->split(m_n_chunks);

    mpi::cout << "    Initialization\n";
    for(auto& part_filename: partition_files) {
        ifstream file;
        if (binary) {
            file.open(part_filename, std::ios::in | std::ios::binary);
        }
        else {
            file.open(part_filename, std::ios::in);
        }

        while(!file.eof()) {
            size_t z, level;
            IndexedObject::index_type idx;

            if (binary) {
                file.read((char*)&idx, sizeof(idx));
                file.read((char*)&level, sizeof(level));
                file.read((char*)&z, sizeof(z));
                file.read((char*)&data_size, sizeof(data_size));
                temp.resize(data_size);
                file.read((char*)temp.data(), data_size * sizeof(double));
            }
            else {
                file >> idx >> level >> z >> data_size;
                temp.resize(data_size);
                for(size_t i = 0; i < data_size; ++i) {
                    file >> temp[i];
                }
            }

            size_t z_base = z / math::pow4(level);
            if (m_base_cells.count(z_base) < 1) {
                continue;
            }

            Cell::Ptr cell = m_base_cells[z_base];

            for(int lev = (int)level - 1; lev >= 0; --lev) {
                size_t z_loc = (z / math::pow4(lev)) % 4;
                if (cell->is_childless()) {
                    cell = nullptr;
                    break;
                }
                cell = cell->get_child_by_local_id(z_loc);
            }

            if (cell != nullptr) {
                cell->data_holder()->deserialize(temp.data());
            }
        }
        file.close();
    }

    mpi::barrier();
    mpi::cout << "  Restart complete\n\n";

    return restart_time;
}