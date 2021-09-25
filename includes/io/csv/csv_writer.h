#pragma once

#include <io/mesh_writer.h>

class NodeList;

class CsvWriter : public MeshWriter {

public:

    explicit CsvWriter(const Configuration &config, double init_time, double t_scale);

    void write_head(Problem* problem) final;

    void write_tail() final;

    void write(Problem* problem, Mesh *mesh, double time) final;

};