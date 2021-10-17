#include <control/configuration.h>
#include <core/face/face.h>
#include <core/cell/cell.h>
#include <core/cell/data_holder.h>
#include <core/mesh/mesh.h>
#include <core/generator/geometry_generator.h>
#include <core/generator/rectangle_generator.h>
#include <problems/surface_problem.h>

SurfaceProblem::SurfaceProblem(const Configuration& config) {

}

uint SurfaceProblem::cell_data_size() const {
    return static_cast<uint>(ceil(double(sizeof(CellData)) / sizeof(double)));
}

CellData SurfaceProblem::get_state(Cell_Ref cell) const {
    CellData res;
    res.vol = cell->data()[0];
    return res;
}

void SurfaceProblem::set_state(Cell_Ref cell, const CellData& value) const {
    cell->data()[0] = value.vol;
}

double is_inside(const Vector3d& r) {
    Vector3d c1 = {0.5, 0.5, 0.0};
    double R1 = 0.3;
    Vector3d c2 = {0.5, 0.0, 0.0};
    double R2 = 0.6;
    if ((r - c1).norm() < R1) {
        if ((r - c2).norm() > R2) {
            return 1.0;
        }
        else {
            return 0.0;
        }
    }
    else {
        return 0.0;
    }
}

void SurfaceProblem::initialization(Mesh *mesh) {
    for (auto cell: *mesh->cells()) {
        CellData data;
        data.vol = is_inside(cell->center());

        set_state(cell, data);
    }
}

double SurfaceProblem::solution_step(Mesh *mesh) {
    return 0.0;
}

AdaptationFlag SurfaceProblem::adaptation_criterion(Cell_Ref cell) {
    double vol = get_state(cell).vol;
    for (auto& face: cell->faces()) {
        auto neib = face->neighbor(cell);

        double neib_vol = get_state(neib).vol;

        if (vol != neib_vol) {
            return AdaptationFlag::SPLIT;
        }
    }

    return AdaptationFlag::COARSE;
}

void SurfaceProblem::split_data(Cell_Ref parent, const vector<Cell_Ptr> &children) {
    auto pdata = get_state(parent);
    for(auto& child: children) {
        child->data_holder()->resize(cell_data_size());
        set_state(child, pdata);
    }
}

void SurfaceProblem::coarse_data(Cell_Ref parent, const vector<Cell_Ptr> &children) {
    auto chdata = get_state(children[0]);
    set_state(parent, chdata);
}

void SurfaceProblem::print_info(const char *tab) const {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << tab << "шташ\n";
}

double SurfaceProblem::get_cell_param(Cell_Ref cell, const string& name) const {
    auto data = get_state(cell);
    if (name == "vol") {
        return data.vol;
    }
    else if (name == "nx") {
        return data.n[0];
    }
    else if (name == "ny") {
        return data.n[1];
    }
    else if (name == "nz") {
        return data.n[2];
    }
    else {
        throw std::runtime_error("Unknown parameter '" + name + "'");
    }
}

double SurfaceProblem::get_integral_param(const std::string& name) const {
    return 0.0;
}
