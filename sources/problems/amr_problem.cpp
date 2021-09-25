#include <control/configuration.h>
#include <core/face/face.h>
#include <core/cell/cell.h>
#include <core/cell/data_holder.h>
#include <core/mesh/mesh.h>
#include <core/generator/geometry_generator.h>
#include <core/generator/rectangle_generator.h>
#include <problems/amr_problem.h>


AmrProblem::AmrProblem(const Configuration& config) {
    m_time = 0.0;
}

uint AmrProblem::cell_data_size() const {
    return 1;
}

void AmrProblem::initialization(Mesh *mesh) {
    mesh->map([this](const NodeList::Part& cells) {
        for(auto cell: cells) {
            set_data(cell, base_function(cell->center(), m_time));
        }
    });
}

void AmrProblem::stage_1(const NodeList::Part& cells) {
    for (auto cell: cells) {

        double u = 0.0;
        double w_sum = 0;
        for(auto& face: cell->faces_list()) {
            auto nei = face->neighbor(cell);

            if (face->flag() != FaceFlag::ORDER) {
                double w_i = 1.0 / (nei->center() - cell->center()).squaredNorm();
                u += w_i * 1.0;
                w_sum += w_i;
            }
            else {
                double w_i = 1.0 / (nei->center() - cell->center()).squaredNorm();
                u += w_i * base_function(nei->center(), m_time);
                w_sum += w_i;
            }
        }

        double sigma = 0.0;
        u = sigma * (u / w_sum) + (1.0 - sigma) * base_function(cell->center(), m_time);

        //double u = base_function(cell->center(), m_time);

        set_data(cell, u);
    }
}

double AmrProblem::solution_step(Mesh *mesh) {
    mesh->sync_data();

    auto calc_time = mesh->map(
            [this](const NodeList::Part &cells) {
                stage_1(cells);
            });

    m_time += 2.0;
    return calc_time;
}

AdaptationFlag AmrProblem::adaptation_criterion(Cell_Ref cell) {
    uint level = (uint)max(0.0, floor(get_data(cell)));

    if (level > cell->level()) {
        return AdaptationFlag::SPLIT;
    }
    else if (level == cell->level()) {
        return AdaptationFlag::NONE;
    }
    else {
        return AdaptationFlag::COARSE;
    }
}

void AmrProblem::split_data(Cell::Ref parent, const vector<Cell::Ptr> &children) {
    for(auto& child: children) {
        child->data_holder()->resize(cell_data_size());
        set_data(child, base_function(child->center(), m_time));
    }
}

void AmrProblem::coarse_data(Cell::Ref parent, const vector<Cell::Ptr> &children) {
    set_data(parent, base_function(parent->center(), m_time));
}

void AmrProblem::print_info(const char *tab) const {
    std::cout << tab << "Time: " << std::setprecision(1) << m_time << "\n";
}

double AmrProblem::get_cell_param(Cell_Ref cell, const string& name) const {
    if (name == "u") {
        return get_data(cell);
    }
    else if (name == "lvl") {
        return cell->level();
    }
    else if (name == "vol") {
        return cell->volume();
    }
    else {
        throw std::runtime_error("Unknown cell parameter '" + name + "'");
    }
}

double AmrProblem::get_integral_param(const std::string& name) const {
    throw std::runtime_error("Has no integral variables");
}

double AmrProblem::get_data(Cell::Ref cell) const {
    return cell->data()[0];
}

void AmrProblem::set_data(Cell::Ref cell, double value) const {
    cell->data()[0] = value;
}

double AmrProblem::base_function(const Vector3d& v0, double t) const {
    double u = 0.0;

    Vector3d v = {fabs(v0[0]), v0[1], v0[2]};

    double omega = 2 * M_PI * t / 1000.0; // omega in [0, 2 pi]

    // Летающий шарик
    {
        uint LEVEL = 5;

        Vector3d C = {500.0, 500.0, 0.0};
        double R = 300.0;
        double r = 70.0;

        double phi = 2 * omega;
        Vector3d c = C;
        c[0] += R * cos(phi);
        c[1] += R * sin(phi);

        if ((v - c).norm() < r) {
            u = LEVEL;
        }
    }

    // Вращающийся крест
    {
        uint LEVEL = 3;

        Vector3d C = {500.0, 500.0, 0.0};
        double H = 25.0;
        double L = 150.0;

        double phi = 3 * omega;

        Vector3d vR = v - C;
        double r = vR.norm();
        double a = atan2(vR[1], vR[0]) + phi;
        vR[0] = r * cos(a);
        vR[1] = r * sin(a);

        if (-H <= vR[0] && vR[0] <= H &&
            -L <= vR[1] && vR[1] <= L) {
            u = LEVEL;
        }
        if (-H <= vR[1] && vR[1] <= H &&
            -L <= vR[0] && vR[0] <= L) {
            u = LEVEL;
        }
    }

    // Понижающие и повышающие полосы
    {
        double H = 100.0;

        double x1 = 300.0;
        double x2 = 700.0;
        if (fabs(v[0] - x1) <= H) {
            u -= 1.0;
        }
        if (fabs(v[0] - x2) <= H) {
            u += 1.0;
        }

        double y1 = 300.0;
        double y2 = 700.0;
        if (fabs(v[1] - y1) <= H) {
            u -= 1.0;
        }
        if (fabs(v[1] - y2) <= H) {
            u += 1.0;
        }
    }


    // Блуждающая горка
    {
        double y = 3 * omega;
        Vector3d c = {1300, 500 * (1.0 + sin(y)), 0.0};
        double s = 300;
        double s2 = s * s;

        double e = 5 * exp(-(v - c).squaredNorm() / s2);
        u += floor(e);
    }

    return u;
}