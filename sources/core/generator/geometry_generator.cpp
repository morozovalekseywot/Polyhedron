//
// Created by 159-mrv on 8/16/18.
//

#include <control/configuration.h>
#include <core/cell/cell.h>
#include <core/generator/structured_decomposition.h>
#include <core/generator/sector_generator.h>
#include <core/generator/sector3d_generator.h>
#include <core/generator/rectangle_generator.h>
#include <core/generator/real_sector_generator.h>
#include <core/generator/cuboid_generator.h>

unique_ptr<GeometryGenerator> GeometryGenerator::create(const Configuration &config) {
    GeometryType type = string_to_geometry(config("mesh", "type").to_string());

    switch (type) {
#if DIM2
        case GeometryType::RECTANGLE:
            return makeUnique<RectangleGenerator>(config);

        case GeometryType::SECTOR:
            return makeUnique<SectorGenerator>(config);

        case GeometryType::REAL_SECTOR:
            return makeUnique<RealSectorGenerator>(config);
#else
        case GeometryType::CUBOID:
            return makeUnique<CuboidGenerator>(config);

        case GeometryType::SECTOR3D:
            return makeUnique<Sector3dGenerator>(config);
#endif
        default:
            throw runtime_error("GeometryGenerator::create() error: Unknown type of geometry");
    }
}

GeometryGenerator::GeometryGenerator(GeometryType type)
        : m_type(type) {

}

GeometryType GeometryGenerator::type() const noexcept {
    return m_type;
}