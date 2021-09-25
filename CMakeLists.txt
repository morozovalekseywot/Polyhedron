cmake_minimum_required(VERSION 3.5)

# =============================== SETTINGS ====================================

project(InfiniMesh)

set(CMAKE_CXX_STANDARD 14)
set(USE_TCMALLOC FALSE)

set(ENABLE_MPI OFF)

# Dimension
add_definitions(-DDIMENSION=2)

#add definition to check some functions in runtime (decrease performance)
#add_definitions(-DFULL_CHECK)


# ============================= FIND PACKAGES =================================

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Threads REQUIRED)
find_package(Boost COMPONENTS system filesystem regex program_options)
if (ENABLE_MPI)
    find_package(HDF5)
    find_package(MPI REQUIRED)
endif()

if (Boost_FOUND)
    message("--   Boost include: ${Boost_INCLUDE_DIR}")
    message("--   Boost libs:    ${Boost_LIBRARIES}")
endif()

if (HDF5_FOUND)
    message(STATUS "HDF5 FOUND")
    message(STATUS "  HDF5 include: ${HDF5_INCLUDE_DIR}")
    message(STATUS "  HDF5 libs:    ${HDF5_LIBRARIES}")
else()
    message(WARNING "HDF5 library is not found, but I can compile project without it.")
endif()



# ================================ COMPILE ====================================

include_directories(includes)
include_directories(libs/json11/include libs/Eigen)
include_directories(${Boost_INCLUDE_DIR})
if (ENABLE_MPI)
    include_directories(${MPI_INCLUDE_PATH})
endif()


# Если тип сборки не указан - использовать релиз
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Добавим флагов компиляции для гарантии
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    # -O3 -- максимальная оптимизация
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -O3")
else()
    # -O0       -- минимальная оптимизация
    # -g        -- что-то связанное с дебагом
    # -rdynamic -- динамическая линковка
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-attributes")

message("-- Build type: ${CMAKE_BUILD_TYPE}")
message("-- C++ compile flags: ${CMAKE_CXX_FLAGS}")


# compile HDFWriter
if (HDF5_FOUND)
    add_definitions(-DUSE_HDF_WRITER)

    include_directories(${HDF5_INCLUDE_DIR})

    set(HDF_WRITER_SOURCES
            includes/io/hdf/hdf_utils.h
            includes/io/hdf/hdf_writer.h
            sources/io/hdf/hdf_writer.cpp
            sources/io/hdf/hdf_utils.cpp)

    add_library(hdfwriter ${HDF_WRITER_SOURCES})
    target_link_libraries(hdfwriter ${HDF5_LIBRARIES})
endif()


set(HEADERS
        includes/allstd.h
        includes/control/configuration.h
        includes/control/calculation_controller.h
        includes/control/mpi_wrapper.h
        includes/control/mpi_commutator.h
        includes/core/indexed_object.h
        includes/core/vertex/vertex.h
        includes/core/face/face.h
        includes/core/face/faces_list.h
        includes/core/cell/data_holder.h
        includes/core/cell/side.h
        includes/core/cell/cell.h
        includes/core/layer/ghost_control.h
        includes/core/layer/border_control.h
        includes/core/layer/layer_control.h
        includes/core/mesh/mesh.h
        includes/core/mesh/balancing_options.h
        includes/core/mesh/partitioner.h
        includes/core/generator/structured_decomposition.h
        includes/core/generator/geometry_generator.h
        includes/core/generator/structured_generator.h
        includes/core/generator/rectangle_generator.h
        includes/core/generator/sector_generator.h
        includes/core/generator/real_sector_generator.h
        includes/core/generator/decomposition.h
        includes/core/generator/orth_decomposition.h
        includes/core/generator/cuboid_generator.h
        includes/core/generator/sector3d_generator.h
        includes/io/console.h
        includes/io/mesh_writer.h
        includes/io/vtk/vtk_writer.h
        includes/io/csv/csv_writer.h
        includes/io/dat/dat_writer.h
        includes/problems/problem.h
        includes/problems/amr_problem.h
        includes/problems/convection_problem.h
        includes/utils/math.h
        includes/utils/geom.h
        includes/utils/interface.h
        includes/utils/stopwatch.h
        includes/utils/memory/kd_tree.h
        includes/utils/memory/node_list.h
        includes/utils/memory/node.h
        includes/utils/memory/list_position.h
        includes/utils/partitioning/splitter.h
        includes/utils/partitioning/base_splitter_1d.h
        includes/utils/partitioning/base_splitter_2d.h
        includes/utils/partitioning/rcb_splitter.h
        includes/utils/partitioning/voronoi_diagram.h)

set(SOURCES
        libs/json11/src/json11.cpp
        sources/control/calculation_controller.cpp
        sources/control/configuration.cpp
        sources/control/mpi_wrapper.cpp
        sources/core/vertex/vertex.cpp
        sources/core/face/face.cpp
        sources/core/face/faces_list.cpp
        sources/core/cell/side.cpp
        sources/core/cell/data_holder.cpp
        sources/core/cell/cell.cpp
        sources/core/layer/ghost_control.cpp
        sources/core/layer/border_control.cpp
        sources/core/layer/layer_control.cpp
        sources/core/mesh/mesh.cpp
        sources/core/mesh/mesh_restart.cpp
        sources/core/mesh/mesh_adaptation.cpp
        sources/core/mesh/mesh_communication.cpp
        sources/core/mesh/mesh_balancing.cpp
        sources/core/mesh/partitioner.cpp
        sources/core/mesh/balancing_options.cpp
        sources/core/generator/structured_decomposition.cpp
        sources/core/generator/geometry_generator.cpp
        sources/core/generator/structured_generator.cpp
        sources/core/generator/sector_generator.cpp
        sources/core/generator/rectangle_generator.cpp
        sources/core/generator/real_sector_generator.cpp
        sources/core/generator/decomposition.cpp
        sources/core/generator/orth_decomposition.cpp
        sources/core/generator/cuboid_generator.cpp
        sources/core/generator/sector3d_generator.cpp
        sources/io/console.cpp
        sources/io/mesh_writer.cpp
        sources/io/vtk/vtk_writer.cpp
        sources/io/csv/csv_writer.cpp
        sources/io/dat/dat_writer.cpp
        sources/problems/problem.cpp
        sources/problems/amr_problem.cpp
        sources/problems/convection_problem.cpp
        sources/utils/math.cpp
        sources/utils/geom.cpp
        sources/utils/interface.cpp
        sources/utils/memory/kd_tree.cpp
        sources/utils/memory/node.cpp
        sources/utils/memory/list_position.cpp
        sources/utils/memory/node_list.cpp
        sources/utils/partitioning/splitter.cpp
        sources/utils/partitioning/base_splitter_1d.cpp
        sources/utils/partitioning/base_splitter_2d.cpp
        sources/utils/partitioning/rcb_splitter.cpp
        sources/utils/partitioning/voronoi_diagram.cpp
        )


add_library(${CMAKE_PROJECT_NAME}_lib STATIC ${SOURCES} ${HEADERS})
target_link_libraries(${CMAKE_PROJECT_NAME}_lib ${Boost_LIBRARIES})
target_link_libraries(${CMAKE_PROJECT_NAME}_lib ${CMAKE_THREAD_LIBS_INIT})
if (ENABLE_MPI)
    target_link_libraries(${CMAKE_PROJECT_NAME}_lib ${MPI_LIBRARIES})
endif()
if (HDF5_FOUND)
    target_link_libraries(${CMAKE_PROJECT_NAME}_lib hdfwriter)
endif()


add_executable(${CMAKE_PROJECT_NAME} main.cpp)
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${boost_INCLUDE_DIRS})
target_link_libraries(${CMAKE_PROJECT_NAME} ${Boost_LIBRARIES})
target_link_libraries(${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_NAME}_lib)


# ============================ BOOST UNIT TESTS ===============================

find_package(Boost COMPONENTS unit_test_framework)
if(Boost_FOUND)

    include_directories(${Boost_INCLUDE_DIR})
    file(GLOB_RECURSE TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} tests/*.cpp)

    foreach(testSrc ${TEST_SRCS})
        message(STATUS "TEST:" ${testSrc})
        get_filename_component(testName ${testSrc} NAME_WE)
        add_executable(${testName} ${testSrc})
        target_link_libraries(${testName} ${Boost_LIBRARIES} ${CMAKE_PROJECT_NAME}_lib)

        set_target_properties(
                ${testName} PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/tests/
        )
        add_test(
                NAME ${testName}
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/tests/
                COMMAND ${PROJECT_SOURCE_DIR}/bin/tests/${testName}
        )
    endforeach()
    enable_testing()
endif()



# ================================= DOXYGEN ===================================

find_package(Doxygen)

if(${DOXYGEN_FOUND})
    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_SOURCE_DIR}/doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc
            COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION doc)
endif()
