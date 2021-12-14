#include <control/configuration.h>
#include <core/face/face.h>
#include <core/cell/cell.h>
#include <core/cell/data_holder.h>
#include <core/mesh/mesh.h>
#include <utils/memory/node_list.h>
#include <core/generator/geometry_generator.h>
#include <core/generator/rectangle_generator.h>
#include <problems/jacobi.h>

Jacobi::Jacobi(const Configuration &config)
{
    figure = surf::Surface("examples/figure/tornado.stl");
    figure.centering();
    figure.rotate({0.0, 0.0, 1.0}, -M_PI_2);
    V0 = 10.0;
}

uint Jacobi::cell_data_size() const
{
    return static_cast<uint>(ceil(double(sizeof(JacobiCellData)) / sizeof(double)));
}

void Jacobi::initialization(Mesh *mesh)
{
    static bool first_time = true;
    if (!first_time) {
        return;
    }
    first_time = false;

    // На стадии инициализации расставляем индексы ячеек, и начальные потенциалы
    // некоторым образом (u1).
    int idx = 0;
    for (auto cell: *mesh->cells())
    {
        JacobiCellData data;
        if (figure.is_inside(cell->center()) || figure.is_inside(cell->center() + Vector3d{0, 0, figure.getMLength() * 0.001}) ||
            figure.is_inside(cell->center() + Vector3d{0, figure.getMLength() * 0.001, 0}) ||
            figure.is_inside(cell->center() + Vector3d{figure.getMLength() * 0.001, 0, 0}))
            data.vol = 1.0;
        else
            data.vol = 0.0;
        data.u1 = -V0 * cell->center()[0];
        data.u2 = -V0 * cell->center()[0];
        data.p = 0.0;
        data.idx = idx++;
        set_state(cell, data);
    }
}

double Jacobi::boundary_function(FaceFlag flag, const Vector3d &vec, const Vector3d &n) const
{
    if (flag == FaceFlag::INFLOW)
        return V0;

    if (flag == FaceFlag::OUTFLOW)
        return -V0;

    return 0.0;
}

void Jacobi::JacobiStage(const NodeList::Part &cells) const
{
    // Во время расчетного шага считаем потенциалы на следующей итерации u2,
    // затем записываем u2 на место u1, считаем градиент u1 методом МНК, то
    // есть получаем компоненты Vx, Vy, Vz
    // mu = S_ab/abs(center_b-center_a)
    // u_a * sum(mu) - sum(mu * u_b) = - граничное условие или 0
    for (auto cell: cells)
    {
        auto zc = get_state(cell);

        // Ячейки тела
        if (get_state(cell).vol > 0.5)
        {
            zc.u1 = 0.0;
            zc.u2 = 0.0;
            zc.p = 0.0;
            zc.v = Vector3d::Zero();

            int count = 0;
            for (auto &face: cell->faces())
            {
                auto neib = face->neighbor(cell);
                auto zn = get_state(neib);

                if (zn.vol > 0.5) {
                    break;
                }

                // Соседняя ячейка в жидкости
                count++;
                zc.u1 += zn.u1;
                zc.u2 += zn.u2;
                zc.v += zn.v;
                zc.p += 0.5 * zn.v.norm() - zn.v[0] * V0;
            }
            if (count > 0)
            {
                zc.u1 /= count;
                zc.u2 /= count;
                zc.p /= count;
                zc.v /= count;

                set_state(cell, zc);
            }
            continue;
        }

        double neib_sum = 0.0, mu_sum = 0.0;
        double cond_sum = 0.0;
        for (auto &face: cell->faces())
        {
            if (face->flag() != FaceFlag::ORDER)
            {
                cond_sum += face->area() * boundary_function(face->flag(), face->center(), -face->normal(cell));
            } else
            {
                auto neib = face->neighbor(cell);
                auto zn = get_state(neib);
                if (zn.vol != zc.vol) {
                    continue;
                }
                double mu = face->area() / (neib->center() - cell->center()).dot(face->normal(neib));
                mu_sum += mu;
                neib_sum += zn.u1 * mu;
            }
        }

        JacobiCellData data = get_state(cell);
        if (mu_sum != 0.0) {
            data.u2 = (cond_sum + neib_sum) / mu_sum;
        }
        else {
            data.u2 = 0.0;
        }

        set_state(cell, data);
    }
}

void Jacobi::VelocityStage(const NodeList::Part &cells) const
{
    // расчёт скорости
    for (auto cell: cells)
    {
        auto cell_c = cell->center();

        JacobiCellData zc = get_state(cell);
        double uc = zc.u1;

        if (get_state(cell).vol > 0.5) {
            zc.v = {0.0, 0.0, 0.0};
            set_state(cell, zc);
            continue;
        }

        Matrix3d A = Matrix3d::Zero();
        Vector3d f = Vector3d::Zero();
        for (auto &face: cell->faces())
        {
            auto face_c = face->center();
            auto normal = -face->normal(cell);

            if (face->flag() != FaceFlag::ORDER)
            {
                auto phi = boundary_function(face->flag(), face_c, normal);

                // Потенциал на грани
                double uf = uc + phi * (cell_c - face_c).norm();

                Vector3d c = face_c - cell_c;
                double w = 1.0 / c.squaredNorm();
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        A(i, j) += w * c(i) * c(j);

                f += w * (uc - uf) * c;
                continue;
            }

            auto zn = get_state(face->neighbor(cell));
            double un = zn.vol != zc.vol ? uc : zn.u1;

            Vector3d c = face->neighbor(cell)->center() - cell->center();
            double w = 1.0 / c.squaredNorm();

            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    A(i, j) += w * c(i) * c(j);

            f += w * (uc - un) * c;
        }

        _2D(A(2, 2) = 1.0;)

        auto B = A.inverse().eval();

        zc.v = B * f;
        set_state(cell, zc);
    }
}

std::array<double, 2> Jacobi::ErrorsStage(const NodeList::Part &cells) const
{
    // расчёт погрешностей
    double eps = 0.0;
    double delta = 0.0;

    for (auto cell: cells)
    {
        auto zc = get_state(cell);

        if (zc.vol > 0.5) {
            zc.eps = 0.0;
            zc.delta = 0.0;
            set_state(cell, zc);
            continue;
        }

        double neib_sum = 0.0, mu_sum = 0.0;
        double cond_sum = 0.0;
        for (auto &face: cell->faces())
        {
            if (face->flag() != FaceFlag::ORDER)
            {
                cond_sum += face->area() * boundary_function(face->flag(), face->center(), -face->normal(cell));
            } else
            {
                auto neib = face->neighbor(cell);
                auto zn = get_state(neib);
                if (zn.vol != zc.vol) {
                    continue;
                }
                double mu = face->area() / (neib->center() - cell->center()).dot(face->normal(neib));
                mu_sum += mu;
                neib_sum += zn.u1 * mu;
            }
        }

        zc.eps = std::abs(zc.u2 - zc.u1);
        zc.delta = std::abs(cond_sum + neib_sum - mu_sum * zc.u1);

        eps = std::max(eps, zc.eps);
        delta = std::max(delta, zc.delta);

        set_state(cell, zc);
    }

    return {eps, delta};
}

void Jacobi::UpdateStage(const NodeList::Part &cells) const
{
    for (auto cell: cells)
    {
        auto data = get_state(cell);
        data.u1 = data.u2;
        data.u2 = 0.0;
        set_state(cell, data);
    }
}

double Jacobi::solution_step(Mesh *mesh)
{
    using std::placeholders::_1;

    mesh->map(std::bind(&Jacobi::JacobiStage, this, _1));

    mesh->map(std::bind(&Jacobi::VelocityStage, this, _1));

    std::vector<double> eps(mesh->n_chunks());
    std::vector<double> delta(mesh->n_chunks());
    mesh->map([this, &eps, &delta](const NodeList::Part &cells)
              {
                  auto res = ErrorsStage(cells);
                  eps[cells.part_id()] = res[0];
                  delta[cells.part_id()] = res[1];
              });
    m_eps = *std::max_element(eps.begin(), eps.end());
    average_eps = 0.0;
    for (auto &e: eps)
        average_eps += e;
    average_eps /= eps.size();

    m_delta = *std::max_element(delta.begin(), delta.end());
    average_delta = 0.0;
    for (auto &d: delta)
        average_delta += d;
    average_delta /= delta.size();

    mesh->map(std::bind(&Jacobi::UpdateStage, this, _1));

    m_time += 1.0;
    return 0;
}

AdaptationFlag Jacobi::adaptation_criterion(const shared_ptr<Cell> &cell)
{
    if (m_time > 0.0)
        return AdaptationFlag::NONE;

    double vol = get_state(cell).vol;
    for (auto &face: cell->faces())
    {
        auto neib = face->neighbor(cell);

        double neib_vol = get_state(neib).vol;

        if (vol != neib_vol)
        {
            return AdaptationFlag::SPLIT;
        }
    }

    return AdaptationFlag::COARSE;
}

void Jacobi::split_data(const shared_ptr<Cell> &parent, const vector<shared_ptr<Cell>> &children)
{
    auto pdata = get_state(parent);
    for (auto &child: children)
    {
        child->data_holder()->resize(cell_data_size());

        auto ch_data = pdata;

        if (figure.is_inside(child->center()) ||
            figure.is_inside(child->center() + Vector3d{0, 0, figure.getMLength() * 0.001}) ||
            figure.is_inside(child->center() + Vector3d{0, figure.getMLength() * 0.001, 0}) ||
            figure.is_inside(child->center() + Vector3d{figure.getMLength() * 0.001, 0, 0}))
            ch_data.vol = 1.0;
        else
            ch_data.vol = 0.0;

        set_state(child, ch_data);
    }
}

void Jacobi::coarse_data(const shared_ptr<Cell> &parent, const vector<shared_ptr<Cell>> &children)
{
    auto p_data = get_state(children[0]);

    if (figure.is_inside(parent->center()) ||
        figure.is_inside(parent->center() + Vector3d{0, 0, figure.getMLength() * 0.001}) ||
        figure.is_inside(parent->center() + Vector3d{0, figure.getMLength() * 0.001, 0}) ||
        figure.is_inside(parent->center() + Vector3d{figure.getMLength() * 0.001, 0, 0}))
        p_data.vol = 1.0;
    else
        p_data.vol = 0.0;

    set_state(parent, p_data);
}

void Jacobi::print_info(const char *tab) const
{
    std::cout << tab << std::scientific << std::setprecision(4)
              << "Eps = " << m_eps << ", "
              << "Average Eps = " << average_eps << ",\n" << tab
              << "Delta = " << m_delta << ", "
              <<"Average Delta = " << average_delta << "\n";
}

double Jacobi::get_cell_param(const shared_ptr<Cell> &cell, const string &name) const
{
    JacobiCellData data = get_state(cell);
    if (name == "u")
        return data.u1;
    else if (name == "vx")
        return data.v.x();
    else if (name == "vy")
        return data.v.y();
    else if (name == "vz")
        return data.v.z();
    else if (name == "vol")
        return data.vol;
    else if (name == "v")
        return data.v.norm();
    else if (name == "p")
        return data.p;
    else if (name == "eps") {
        return data.eps;
    }
    else if (name == "delta") {
        return data.delta;
    }
    else
        throw std::runtime_error("Unknown parameter '" + name + "'");
}

double Jacobi::get_integral_param(const string &name) const
{
    if (name == "eps")
    {
        return m_eps;
    } else if (name == "delta")
    {
        return m_delta;
    } else
    {
        throw std::runtime_error("Unknown integral parameter '" + name + "'");
    }
}

void Jacobi::set_state(const shared_ptr<Cell> &cell, const JacobiCellData &value) const
{
    memcpy(cell->data(), &value, sizeof(JacobiCellData));
}

JacobiCellData Jacobi::get_state(const shared_ptr<Cell> &cell) const
{
    JacobiCellData res;
    memcpy(&res, cell->data(), sizeof(JacobiCellData));
    return res;
}
