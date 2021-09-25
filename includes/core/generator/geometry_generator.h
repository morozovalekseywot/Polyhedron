//
// Created by 159-mrv on 3/28/18.
//

#ifndef INFINIMESH_GEOMETRY_GENERATOR_H
#define INFINIMESH_GEOMETRY_GENERATOR_H

#include <allstd.h>
#include <core/cell/side.h>

class Cell;
class Face;
class Vertex;
class Configuration;
class Decomposition;

enum class Side;

enum class GeometryType : int {
//#if DIM2
    SECTOR, RECTANGLE, REAL_SECTOR,
#if DIM3
    CUBOID, SECTOR3D
#endif
};

inline GeometryType string_to_geometry(const string& str) {
#if DIM2
    if(str == "sector") {
        return GeometryType::SECTOR;
    }
    else if(str == "rectangle") {
        return GeometryType::RECTANGLE;
    }
    else if (str == "real_sector" ){
        return GeometryType::REAL_SECTOR;
    } else {
        throw std::runtime_error("Error! No such geometry (" + str + ")!");
    }
#else
    if (str == "cuboid") {
        return GeometryType::CUBOID;
    }
    else if (str == "sector3d") {
        return GeometryType::SECTOR3D;
    }
    else {
        throw std::runtime_error("Error! No such geometry (" + str + ")!");
    }
#endif
}

class GeometryGenerator {

protected:
    using Vertex_Ptr = shared_ptr<Vertex>;
    using Face_Ptr = shared_ptr<Face>;
    using Face_Ref = const shared_ptr<Face>&;
    using Cell_Ptr = shared_ptr<Cell>;
    using Cell_Ref = const shared_ptr<Cell>&;

public:
    using UPtr = unique_ptr<GeometryGenerator>;

    static GeometryGenerator::UPtr create(const Configuration &config);

    explicit GeometryGenerator(GeometryType type);

    GeometryType type() const noexcept;

    virtual vector<Cell_Ptr> create_cells(Decomposition *decomp) const = 0;

    virtual Cell_Ptr create_base_cell(size_t z) const = 0;

    virtual Face_Ptr next_face_by_side(Cell_Ref cell, Side side) const = 0;

    virtual size_t neighbor_z(Cell_Ref cell, Side side) const = 0;

    virtual size_t n_cells() const = 0;

    virtual void print_info(const string &tab) const = 0;

protected:
    GeometryType m_type;
};

#endif //INFINIMESH_GEOMETRY_GENERATOR_H
