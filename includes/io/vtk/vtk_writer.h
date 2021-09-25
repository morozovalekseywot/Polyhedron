//
// Created by 159-egi on 2/17/17.
//

#ifndef INFINIMESH_VTK_WRITER_H
#define INFINIMESH_VTK_WRITER_H

#include <io/mesh_writer.h>

class NodeList;

class VtkWriter : public MeshWriter {

public:
    using SPtr = shared_ptr<VtkWriter>;

    explicit VtkWriter(const Configuration &config, double init_time, double t_scale);

    void write_head(Problem* problem) override;
    void write_tail() override;
    void write(Problem* problem, Mesh *mesh, double time) override;

    static void write(const NodeList &cells, const string &filename, Problem* problem = nullptr,
                      const std::vector<std::string>& params = {});

    static const std::vector<std::string>& default_params();

private:
    static void write_binary_xml(Problem* problem, const NodeList &cells, const string &filename,
                                 const std::vector<std::string>& params);
    static void write_binary_header(Problem* problem, const NodeList &cells, std::ofstream &os,
                                    const std::vector<std::string>& params);
    static void write_binary_data(Problem* problem, const NodeList &cells, std::ofstream &os,
                                  const std::vector<std::string>& params);
    static void write_binary_footer(std::ofstream& os);

    void write_binary_xml(Problem* problem, const NodeList &cells);
    void write_pvd_step(double time);

    static std::vector<std::string> m_default_variables;
};

#endif //INFINIMESH_VTK_WRITER_H
