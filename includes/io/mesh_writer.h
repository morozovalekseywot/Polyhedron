//
// Created by 159-mrv on 10/26/17.
//

#ifndef AMR_MESH_WRITER_H
#define AMR_MESH_WRITER_H

#include <allstd.h>

class Mesh;
class Problem;

class MeshWriter {

public:
    explicit MeshWriter(const Configuration &config, double init_time, double t_scale = 1.0);

    static unique_ptr<MeshWriter> create(const Configuration &config, double init_time, double t_scale = 1.0);

    virtual void write_head(Problem* problem) = 0;
    virtual void write_tail() = 0;
    virtual void write(Problem* problem, Mesh *mesh, double time) = 0;

    static std::string pt(int rank = -1);

protected:

    bool time_to_write(double time);

    string m_filename;
    string m_fullname;
    string m_tempname;

    vector<double> m_timecodes;
    vector<double> m_frequencies;
    double m_next_write;
    double m_time_delta;
    size_t m_step;

    std::vector<std::string> m_variables;
};

#endif //AMR_MESH_WRITER_H
