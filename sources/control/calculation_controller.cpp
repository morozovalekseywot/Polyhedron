#include <control/calculation_controller.h>
#include <control/configuration.h>
#include <control/mpi_wrapper.h>

#include <core/indexed_object.h>
#include <core/face/faces_list.h>
#include <core/cell/cell.h>
#include <core/mesh/mesh.h>

#include <io/mesh_writer.h>
#include <io/vtk/vtk_writer.h>

#include <problems/problem.h>
#include <core/cell/side.h>
#include <utils/stopwatch.h>

using namespace std::chrono;

CalculationController::CalculationController(const string& config_filename, const string& checkpoint_filename) {
    Configuration config(config_filename);

    // Ограничители расчета
    m_max_step = std::numeric_limits<size_t>::max();
    m_max_time = std::numeric_limits<double>::max();
    if (config.exist("calculation", "max_time")) {
        m_max_time = config("calculation", "max_time").to_double();
    }
    if (config.exist("calculation", "max_step")) {
        m_max_step = config("calculation", "max_step").to_uint();
    }

    if (config.exist("calculation", "command_file")) {
        m_command_filename = config("calculation", "command_file").to_string();
    }

    // Адаптация
    m_adaptation = false;
    if (config.exist("mesh", "adaptation")) {
        m_adaptation = config("mesh", "adaptation").to_bool();
    }

    m_adapt_frequency = 1;
    if (config.exist("mesh", "adaptation_criteria", "frequency")) {
        m_adapt_frequency = config("mesh", "adaptation_criteria", "frequency").to_uint();
    }

    m_balancing_frequency = 50;
    if (config.exist("mesh", "balancing", "frequency")) {
        m_balancing_frequency = config("mesh", "balancing", "frequency").to_uint();
    }

    // Список чекпоинтов
    m_checkpoint_prefix = "checkpoint_filename";
    if (config.exist("io", "checkpoint_prefix"))
        m_checkpoint_prefix = config("io", "checkpoint_prefix").to_string();
    if (config.exist("io", "output_directory"))
        m_checkpoint_prefix = config("io", "output_directory").to_string() + "/" + m_checkpoint_prefix;

    if (config.exist("io", "checkpoints")) {
        for (auto &checkpoint: config("io", "checkpoints").to_array()) {
            m_checkpoints.push_back(checkpoint.to_double());
        }
    }
    m_checkpoints.push_back(std::numeric_limits<double>::max());
    m_checkpoint_binary = true;
    if (config.exist("io", "checkpoint_format")) {
        if (config("io", "checkpoint_format").to_string() == "binary") {
            m_checkpoint_binary = true;
        }
        else if (config("io", "checkpoint_format").to_string() == "ascii") {
            m_checkpoint_binary = false;
        }
        else {
            throw runtime_error("Error: Unknown checkpoint file format. Available: 'binary'/'ascii'");
        }
    }

    m_problem = Problem::create(config);

    double start_time = 0.0;
    if (checkpoint_filename.empty()) {
        m_mesh = Mesh::create(config);
        m_mesh->extend_cell_data(m_problem->cell_data_size());

        if (m_adaptation) {
            mpi::cout << "  Initial adaptation started\n";
            mpi::cout << "    Refine to\n";
            for (size_t i = 0; i <= m_mesh->max_level(); ++i) {
                mpi::barrier();
                mpi::cout << "      level " << i + 1 << "\n";
                m_problem->initialization(m_mesh.get());
                m_mesh->sync_data();
                m_mesh->adaptation(m_problem.get());
                //m_mesh->dbg_write(m_problem.get());
                mpi::barrier();
                m_mesh->load_balance();
                //m_mesh->dbg_write(m_problem.get());
            }
            mpi::cout << "  Initial adaptation complete\n\n";
        }

        mpi::cout << "  Set initial conditions\n";
        m_problem->initialization(m_mesh.get());

        mpi::barrier();
        mpi::cout << "  OK. RUN!\n\n";
        std::cout.flush();

        start_time = 0.0;
    } else {
        m_mesh = Mesh::create(config, checkpoint_filename);
        //m_mesh->dbg_write(m_problem.get());

        std::ifstream checkpoint_file(checkpoint_filename, std::ios::in);
        char word[6];
        checkpoint_file >> word >> start_time;
        checkpoint_file.close();
    }

    m_problem->set_time(start_time);

    if (config.exist("io", "writers")) {
        if (config("io", "writers").is_array()) {
            for (auto& io_config: config("io", "writers").to_array()) {
                m_writers.emplace_back(MeshWriter::create(io_config, start_time));
            }
        }
        else {
            throw std::runtime_error("Subsection 'io:writers' is not an array");
        }
    }
    else {
        m_writers.emplace_back(MeshWriter::create(config("io"), start_time));
    }
}

CalculationController::~CalculationController() = default;

string extended_time_format(long seconds) {
    long minutes = seconds / 60   % 60;
    long hours   = seconds / 3600 % 24;
    long days    = seconds / 86400;

    stringstream ss;
    if (days > 0) {
        ss << std::setw(2) << days << " d ";
    }
    else {
        ss << "     ";
    }
    if (hours > 0 || days > 0) {
        ss << std::setw(2) << hours << " h ";
    }
    else {
        ss << "     ";
    }
    ss << std::setw(2) << minutes << " m ";
    ss << std::setw(2) << seconds % 60 << " s";

    return ss.str();
}

void CalculationController::print_step_info(double step_time, long full_time) {
    if(mpi::is_master()) {
        std::cout << " +======================================================================+\n";
        std::cout << " |                           Calculation info                           |\n";
        std::cout << " +----------------------------------------------------------------------+\n";

        std::cout << std::fixed;
        std::cout << "  MAIN LOOP ITERATION: " << std::setw(7) << m_step << "                 "
                  << "  CPU SECONDS: " << std::setw(10) << std::setprecision(5) << step_time << "\n";
        std::cout << "                                               "
                  << "  FULL TIME:   " << std::setw(10) << full_time << "\n";
        std::cout << "                                                     "
                  << extended_time_format(full_time) << "\n";

        std::cout << " ------------------------------------------------------------------------\n";

        m_problem->print_info("  ");

        std::cout << " ------------------------------------------------------------------------\n";
        std::cout.flush();
    }
}

void CalculationController::print_time_stat(long calc_time,
                                      long write_time,
                                      long solve_time,
                                      long adapt_time,
                                      long balance_time) {
    if (mpi::is_master()) {
        std::cout << "\n";
        std::cout.flush();

        std::cout << "  Calculation statistics.\n";

        std::cout << "  Full calculation time:  "
                  << extended_time_format(calc_time / 1000)
                  << "  (" << std::setw(10) << calc_time << " ms)\n";

        std::cout << "    Solve time:           "
                  << extended_time_format(solve_time / 1000)
                  << "  (" << std::setw(10) << solve_time << " ms"
                  << std::fixed << setprecision(1) << std::setw(8)
                  << 100.0 * solve_time / calc_time << " %)\n";

        std::cout << "    Write time:           "
                  << extended_time_format(write_time / 1000)
                  << "  (" << std::setw(10) << write_time << " ms"
                  << std::fixed << setprecision(1) << std::setw(8)
                  << 100.0 * write_time / calc_time << " %)\n";

        std::cout << "    Adaptation time:      "
                  << extended_time_format(adapt_time / 1000)
                  << "  (" << std::setw(10) << adapt_time << " ms"
                  << std::fixed << setprecision(1) << std::setw(8)
                  << 100.0 * adapt_time / calc_time << " %)\n";

        std::cout << "    Balancing time:       "
                  << extended_time_format(balance_time / 1000)
                  << "  (" << std::setw(10) << balance_time << " ms"
                  << std::fixed << setprecision(1) << std::setw(8)
                  << 100.0 * balance_time / calc_time << " %)\n";

        std::cout.flush();
    }
}

CalculationController::CommandFlag CalculationController::exec_command() {
    if (m_command_filename.empty())
        return CommandFlag::NONE;

    ifstream command_file(m_command_filename, std::ios::in);
    if (!command_file) {
        return CommandFlag::NONE;
    }

    string command;
    getline(command_file, command);

    std::cout << "FIND COMMAND. IT IS |" << command << "|\n";

    if (command == "exit") {
        command_file.close();
        return CommandFlag::EXIT;
    }
    if (command == "write") {
        string prefix;
        getline(command_file, prefix);
        string filename = prefix + MeshWriter::pt() + ".vtu";
        VtkWriter::write(*m_mesh->cells(), filename, m_problem.get());
        command_file.close();
        mpi::barrier();

        if (mpi::is_master()) {
            ofstream answer_file(m_command_filename, std::ios::out | std::ios::trunc);
            answer_file << "answer\nsuccess write '" << prefix << "'";
            answer_file.close();
        }

        return CommandFlag::NONE;
    }

    command_file.close();
    return CommandFlag::NONE;
}

void CalculationController::run() {
    m_step = 1;

    // Секундомеры
    Stopwatch sw_main;
    Stopwatch sw_write;
    Stopwatch sw_solve;
    Stopwatch sw_adapt;
    Stopwatch sw_balance;

    size_t chkp_n = 0;
    // Это если мы с чекпоинта стартанули
    double current_time = m_problem->current_time();
    if (current_time > 0) {
        while (m_checkpoints[chkp_n] < current_time) {
            ++chkp_n;
        }
    }

    for(auto& writer: m_writers) {
        writer->write_head(m_problem.get());
    }

    double integral_proc_time = 0.0;

    while (current_time < m_max_time && m_step <= m_max_step) {
        sw_main.start();

        // Выполнить runtime комманду
        if (m_step % 10 == 0) {
            auto flag = exec_command();
            if (flag == CommandFlag::EXIT) {
                std::cout << "  Canceled due to runtime command 'exit'\n";
                break;
            }
        }


        // Сохраним сетку
        sw_write.start();
        for(auto& writer: m_writers) {
            writer->write(m_problem.get(), m_mesh.get(), current_time);
        }

        // Сделаем checkpoint
        if (current_time >= m_checkpoints[chkp_n]) {
            string path = m_checkpoint_prefix + "_" + static_cast<char>(chkp_n + 'a');
            m_mesh->checkpoint(path, current_time, m_checkpoint_binary);
            chkp_n = std::min(chkp_n + 1, m_checkpoints.size() - 1);
        }
        sw_write.stop();


        // Шаг решения
        sw_solve.start();
        auto proc_time = solve();
        //perf_test();
        sw_solve.stop();


        // Адаптация
        sw_adapt.start();
        if (m_step % m_adapt_frequency == 0) {
            m_mesh->adaptation(m_problem.get());
        }
        sw_adapt.stop();


        // Балансировка нагрузки
        sw_balance.start();
        integral_proc_time += proc_time;
        if (m_step % m_balancing_frequency == 0) {
            //m_mesh->load_balance(m_mesh->cells()->size());
            m_mesh->load_balance(0.0) ;// integral_proc_time);
            integral_proc_time = 0.0;
            //m_mesh->dbg_write(m_problem.get());
        }
        sw_balance.stop();


        double step_time = sw_main.stop();

        // Вывести информацию
        print_step_info(step_time, (long)sw_main.seconds());
        if (m_step % 1000 == 0) {
            print_time_stat(sw_main.milliseconds_l(),
                            sw_write.milliseconds_l(),
                            sw_solve.milliseconds_l(),
                            sw_adapt.milliseconds_l(),
                            sw_balance.milliseconds_l()
            );
        }

        ++m_step;
        current_time = m_problem->current_time();
    }

    // Запишем последнее состояние
    for(auto& writer: m_writers) {
        writer->write(m_problem.get(), m_mesh.get(), current_time);
        writer->write_tail();
    }

    sw_main.stop();

    print_time_stat(sw_main.milliseconds_l(),
                    sw_write.milliseconds_l(),
                    sw_solve.milliseconds_l(),
                    sw_adapt.milliseconds_l(),
                    sw_balance.milliseconds_l()
    );
}

double CalculationController::solve() {
    auto proc_time = m_problem->solution_step(m_mesh.get());
    return proc_time;
}

void CalculationController::perf_test() {
    static const uint loops = 30;

    vector<uint> thread_nums = {1u, 2u, 4u, 8u, 16u, 24u, 32u};
    vector<double> times_mean(thread_nums.size());
    vector<double> times_disp(thread_nums.size());
    vector<double> ratio(thread_nums.size());
    vector<double> ratio_disp(thread_nums.size());

    for(uint i = 0; i < thread_nums.size(); ++i) {
        m_mesh->cells()->split(thread_nums[i]);
        std::cout << "  start on " << thread_nums[i] << " CPUs\n";

        vector<double> times(loops, 0.0);
        for (uint j = 0; j < loops; ++j) {
            times[j] = m_problem->solution_step(m_mesh.get());
        }

        times_mean[i] = math::mean(times);
        times_disp[i] = math::std_deviation(times, times_mean[i]);

        ratio[i] = times_mean[0] / times_mean[i];
        ratio_disp[i] = ratio[i] * sqrt(pow(times_disp[0] / times_mean[0], 2) +
                                        pow(times_disp[i] / times_mean[i], 2));

        std::cout << "    avg_time:    " << std::fixed << std::setprecision(1)
                  << 1000.0 * times_mean[i] << " ± " << 1000.0 * times_disp[i] << " ms\n";
        std::cout << "    performance: " << std::fixed << std::setprecision(1)
                  << ratio[i] << " ± " << ratio_disp[i] << "\n";
    }
}
