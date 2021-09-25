//
// Created by 159-mrv on 2/15/19.
//

#include <control/configuration.h>
#include <core/generator/decomposition.h>
#include <core/generator/geometry_generator.h>
#include <core/generator/sector_generator.h>
#include <core/generator/real_sector_generator.h>
#include <core/generator/rectangle_generator.h>
#include <core/generator/orth_decomposition.h>

unique_ptr<Decomposition> Decomposition::create(const Configuration& config,
                                                GeometryGenerator* generator) {
#if DIM2
    if (generator->type() == GeometryType::SECTOR) {
        return makeUnique<OrthDecomposition>(config, generator);
    } else if (generator->type() == GeometryType::REAL_SECTOR) {
        return makeUnique<OrthDecomposition>(config, generator);
    } else if (generator->type() == GeometryType::RECTANGLE) {
        return makeUnique<OrthDecomposition>(config, generator);
    } else {
        throw runtime_error("Decomposition error! Unknown type of geometry");
    }
#else
    if (generator->type() == GeometryType::CUBOID) {
        return makeUnique<OrthDecomposition>(config, generator);
    } else if (generator->type() == GeometryType::SECTOR3D) {
        return makeUnique<OrthDecomposition>(config, generator);
    } else {
        throw runtime_error("Decomposition error! Unknown type of geometry");
    }
#endif
}