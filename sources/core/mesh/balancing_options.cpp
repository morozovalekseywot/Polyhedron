//
// Created by 159-mrv on 10/16/18.
//

#include <control/configuration.h>
#include <control/mpi_wrapper.h>
#include <core/mesh/balancing_options.h>

BalancingOptions::BalancingOptions(const Configuration &config) {
    if (!config.exist("mesh", "balancing") || mpi::size() < 2) {
        type = BalancingType::NONE;
        cells = CellsSet::BASE;
        proc_per_x = -1;
        proc_per_y = -1;
        return ;
    }

    Configuration options = config("mesh", "balancing");

    string type_s = options("type").to_string();
    if (type_s == "none") {
        type = BalancingType::NONE;
        proc_per_x = -1;
        proc_per_y = -1;

    } else if (type_s == "X") {
        type = BalancingType::X;
        proc_per_y = 1;
        if (options.exist("proc_per_y")) {
            proc_per_y = options("proc_per_y").to_int();
            if (mpi::size() % proc_per_y != 0) {
                proc_per_y = 1;
                if (mpi::is_master()) {
                    std::cerr << "  Warning: Configuration file includes parameters 'mesh.balancing.proc_per_y',\n";
                    std::cerr << "           but mpi::size() is not a multiple of the value.\n";
                    std::cerr << "           The value will be ignored.\n\n";
                }
            }
        }
        proc_per_x = mpi::size() / proc_per_y;

    } else if (type_s == "Y") {
        type = BalancingType::Y;
        proc_per_x = 1;
        if (options.exist("proc_per_x")) {
            proc_per_x = options("proc_per_x").to_int();
            if (mpi::size() % proc_per_x != 0) {
                proc_per_x = 1;
                if (mpi::is_master()) {
                    std::cerr << "  Warning: Configuration file includes parameters 'mesh.balancing.proc_per_x',\n";
                    std::cerr << "           but mpi::size() is not a multiple of the value.\n";
                    std::cerr << "           The value will be ignored.\n\n";
                }
            }
        }
        proc_per_y = mpi::size() / proc_per_x;

    } else if (type_s == "XY") {
        type = BalancingType::XY;
        if (options.exist("proc_per_x") || options.exist("proc_per_y")) {
            if (mpi::is_master()) {
                std::cerr << "  Warning: Configuration file includes parameter 'mesh.balancing.proc_per_x' or/and 'mesh.balancing.proc_per_y'\n";
                std::cerr << "           which has no sense for two dimensional decomposition.\n";
                std::cerr << "           The values will be ignored.\n\n";
            }
        }
        proc_per_x = -1;
        proc_per_y = -1;


    } else {
        throw runtime_error(R"(Error: Incorrect value of mesh.balancing.type, avaliable: "none", "X", "Y", "XY".)");
    }

    string cells_s = options("cells").to_string();
    if (options("cells").to_string() == "base") {
        cells = CellsSet::BASE;
    } else if (options("cells").to_string() == "leaf") {
        cells = CellsSet::LEAF;
    } else {
        throw runtime_error(R"(Error: Incorrect value of mesh.balancing.cells, avaliable: "base", "leaf".)");
    }
}