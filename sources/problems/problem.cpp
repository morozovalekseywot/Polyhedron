#include <control/configuration.h>
#include <problems/problem.h>
#include <problems/amr_problem.h>
#include <problems/convection_problem.h>
#include <problems/game_of_life.h>
#include <problems/surface_problem.h>


Problem::Problem()
    : m_time(0.0) {
}

unique_ptr<Problem> Problem::create(const Configuration &config) {
    if (config.exist("problem", "name")) {
        auto problem = config("problem", "name").to_string();
        if (problem == "AMR") {
            return makeUnique<AmrProblem>(config);
        }
        else if (problem == "Convection") {
            return makeUnique<ConvectionProblem>(config);
        }
        else if (problem == "GameOfLife") {
            return makeUnique<GameOfLife>(config);
        }
        else if (problem == "Surface") {
            return makeUnique<SurfaceProblem>(config);
        }
        else {
            throw runtime_error("Unknown problem '" + problem + "' in config file.");
        }
    }
    else {
        throw runtime_error("Problem is not defined.");
    }
}

double Problem::current_time() const {
    return m_time;
}

void Problem::set_time(double time) {
    m_time = time;
}
