#include <control/configuration.h>
#include <core/face/face.h>
#include <core/cell/cell.h>
#include <core/cell/data_holder.h>
#include <core/mesh/mesh.h>
#include <core/generator/geometry_generator.h>
#include <core/generator/rectangle_generator.h>
#include <problems/game_of_life.h>

GameOfLife::GameOfLife(const Configuration& config) {

}

uint GameOfLife::cell_data_size() const {
    return 2;
}

double GameOfLife::get_old_state(Cell_Ref cell) const {
    return cell->data()[0];
}

void GameOfLife::set_old_state(Cell_Ref cell, double value) const {
    cell->data()[0] = value;
}

double GameOfLife::get_new_state(Cell_Ref cell) const {
    return cell->data()[1];
}

void GameOfLife::set_new_state(Cell_Ref cell, double value) const {
    cell->data()[1] = value;
}

void GameOfLife::initialization(Mesh *mesh) {
    srand(time(0));

    for(auto cell: *mesh->cells()) {
        double val = double(rand() % 2);
        set_old_state(cell, val);

        set_new_state(cell, 0.0);
    }
}

double GameOfLife::solution_step(Mesh *mesh) {
    for (auto cell: *mesh->cells()) {
        double around = 0.0;

        for (auto& face: cell->faces()) {
            auto neib = face->neighbor(cell);

            around += get_old_state(neib);
        }

        double val = get_old_state(cell);

        if (val == 0.0) {
            if (around == 3.0) {
                set_new_state(cell, 1.0);
            }
            else {
                set_new_state(cell, 0.0);
            }
        }
        else {
            if (2.0 <= around && around <= 3.0) {
                set_new_state(cell, 1.0);
            } else {
                set_new_state(cell, 0.0);
            }
        }
    }

    for (auto cell: *mesh->cells()) {
        double val = get_new_state(cell);

        set_old_state(cell, val);
        set_new_state(cell, 0.0);
    }

    m_sum = 0.0;
    for (auto cell: *mesh->cells()) {
        m_sum += get_old_state(cell);
    }
    m_sum /= mesh->cells()->size();

    m_time += 1.0;
    return 0.0;
}

AdaptationFlag GameOfLife::adaptation_criterion(Cell_Ref cell) {
    return AdaptationFlag::NONE;
}

void GameOfLife::print_info(const char *tab) const {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << tab << "Осталось " << 100.0*m_sum << "% выживших\n";
}

double GameOfLife::get_cell_param(Cell_Ref cell, const string& name) const {
    if (name == "u") {
        return get_old_state(cell);
    }
    else if (name == "v") {
        return get_new_state(cell);
    }
    else {
        throw std::runtime_error("Unknown parameter '" + name + "'");
    }
}

double GameOfLife::get_integral_param(const std::string& name) const {
    return 0.0;
}
