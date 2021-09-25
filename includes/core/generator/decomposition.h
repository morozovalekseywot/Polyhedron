//
// Created by 159-mrv on 2/15/19.
//

#ifndef INFINIMESH_DECOMPOSITION_H
#define INFINIMESH_DECOMPOSITION_H

#include <allstd.h>

class Configuration;
class Vertex;
class Cell;
class GeometryGenerator;


class Decomposition {
public:

    static unique_ptr<Decomposition> create(const Configuration& config, GeometryGenerator* geomerty);

    virtual int rank(const Vector3d& vertex) const = 0;

    virtual void set_workload(double workload) = 0;

    virtual void sync_data() = 0;

    virtual void update() = 0;

    virtual ~Decomposition() = default;

};

#endif //INFINIMESH_DECOMPOSITION_H
