#include <string>
#include <boost/program_options.hpp>

#include <control/calculation_controller.h>
#include <control/mpi_wrapper.h>

namespace po = boost::program_options;

po::variables_map parse_options(int argc, char ** argv) {
    po::options_description description("InfiniMesh");

    description.add_options()
            ("config,c", po::value<std::string>()->default_value("../examples/life.json"),
             "configuration file path (default: ../examples/demo.json")
            ("restart,r", po::value<std::string>()->default_value(""),
             "path to checkpoint file");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(description).allow_unregistered().run(), vm);
    po::notify(vm);

    return vm;
}

int main(int argc, char** argv) {
    // Инициализация mpi
    // Обязательно в начале
    mpi::init(argc, argv);

    // Парсим командную строку
    auto vm = parse_options(argc, argv);
    std::string config_filename = vm["config"].as<std::string>();
    std::string checkpoint_filename = vm["restart"].as<std::string>();

    // Создание и запуск CalculationController
    CalculationController controller(config_filename, checkpoint_filename);
    controller.run();

    mpi::finilize();
    return 0;
}
