//
// Created by 159-mrv on 10/26/17.
//

#ifndef AMR_HDF_WRITER_H
#define AMR_HDF_WRITER_H

#include <io/mesh_writer.h>

class HdfWriter : public MeshWriter {

public:
    using SPtr = shared_ptr<HdfWriter>;

    explicit HdfWriter(const Configuration &config, double init_time, double t_scale);

    void write_head(Problem* problem) override;
    void write_tail() override;
    void write(Problem* problem, Mesh *mesh, double time) override;
};


#endif //AMR_HDF_WRITER_H
