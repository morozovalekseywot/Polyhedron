#include <control/configuration.h>
#include <core/face/face.h>
#include <core/cell/cell.h>
#include <core/cell/data_holder.h>
#include <core/mesh/mesh.h>
#include <core/generator/geometry_generator.h>
#include <core/generator/rectangle_generator.h>
#include <problems/convection_problem.h>


ConvectionProblem::ConvectionProblem(const Configuration& config) {
    m_time = 0.0;
}

uint ConvectionProblem::cell_data_size() const {
    return 2;
}

void ConvectionProblem::initialization(Mesh *mesh) {
    Vector3d c = {0.5, 1.5, 0.0};
    double R = 0.3;

    for (auto cell: mesh->all_chunks()) {
        double u = 0.0;
        if ((cell->center() - c).norm() < R) {
            u = 1.0;
        }
        set_old_state(cell, u);
    }
}

void ConvectionProblem::compute_dt(const NodeList::Part& cells) {
    m_dt = std::numeric_limits<double>::max();
    for (auto cell: cells) {
        double dx = std::sqrt(cell->volume());
        double V = velocity(cell->center(), m_time).norm();
        double dt = 0.5 * dx / V;

        m_dt = std::min(m_dt, dt);
    }
}

void ConvectionProblem::make_step(const NodeList::Part& cells) {
    for (auto cell: cells) {
        double zc = get_old_state(cell);

        double fluxes = 0.0;
        for (auto& face: cell->faces()) {
            if (face->flag() != FaceFlag::ORDER) {
                continue;
            }

            auto neib = face->neighbor(cell);
            double zn = get_old_state(neib);

            Vector3d normal = face->normal(neib);
            Vector3d face_c = face->center();
            Vector2d V = velocity(face_c, m_time);

            double Vn = V[0]*normal[0] + V[1]*normal[1];
            double Vnp = std::max(0.0, Vn);
            double Vnm = std::min(0.0, Vn);

            fluxes += (Vnp * zc + Vnm * zn) * face->area();
        }

        double z_new = zc - m_dt * fluxes / cell->volume();

        set_new_state(cell, z_new);
    }
}

void ConvectionProblem::update(const NodeList::Part& cells) {
    for(auto cell: cells) {
        set_old_state(cell, get_new_state(cell));
    }
}

double ConvectionProblem::solution_step(Mesh *mesh) {
    mesh->sync_data();

    compute_dt(mesh->all_chunks());

    make_step(mesh->all_chunks());

    update(mesh->all_chunks());

    m_time += m_dt;

    return 0.0;
}

AdaptationFlag ConvectionProblem::adaptation_criterion(Cell::Ref cell) {
    double u = get_old_state(cell);

    if (u > 0.1) {
        return AdaptationFlag::SPLIT;
    }
    else if (u > 0.05) {
        return AdaptationFlag::NONE;
    }
    else {
        return AdaptationFlag::COARSE;
    }
}

void ConvectionProblem::split_data(Cell::Ref parent, const vector<Cell::Ptr> &children) {
    for(auto& child: children) {
        child->data_holder()->resize(cell_data_size());
        set_old_state(child, get_old_state(parent));
        set_new_state(child, get_new_state(parent));
    }
}

void ConvectionProblem::coarse_data(Cell::Ref parent, const vector<Cell::Ptr> &children) {
    double u1 = 0.0;
    double u2 = 0.0;
    for (auto& child: children) {
        u1 += get_old_state(child);
        u2 += get_new_state(child);
    }
    set_old_state(parent, 0.25*u1);
    set_new_state(parent, 0.25*u2);
}

void ConvectionProblem::print_info(const char *tab) const {
    std::cout << tab << "Time: " << std::setprecision(6) << m_time
              << "\t\tdt: " << m_dt << "\n";
}

double ConvectionProblem::get_cell_param(Cell_Ref cell, const string& name) const {
    if (name == "u") {
        return get_old_state(cell);
    }
    else if (name == "vx") {
        return velocity(cell->center(), m_time)[0];
    }
    else if (name == "vy") {
        return velocity(cell->center(), m_time)[1];
    } else if (name == "lvl") {
        return cell->level();
    }
    else {
        throw std::runtime_error("Unknown cell parameter '" + name + "'");
    }
}

double ConvectionProblem::get_integral_param(const std::string& name) const {
    throw std::runtime_error("Has no integral variables");
}

double ConvectionProblem::get_old_state(Cell::Ref cell) const {
    return cell->data()[0];
}

double ConvectionProblem::get_new_state(Cell::Ref cell) const {
    return cell->data()[1];
}

void ConvectionProblem::set_old_state(Cell::Ref cell, double value) const {
    cell->data()[0] = value;
}

void ConvectionProblem::set_new_state(Cell::Ref cell, double value) const {
    cell->data()[1] = value;
}

Vector2d ConvectionProblem::velocity(const Vector3d& r, double t) const {
    double vx = 1.0;
    //double vy = 0.0;
    double vy = std::sin(3.0*r[0]);
    return {vx, vy};
}