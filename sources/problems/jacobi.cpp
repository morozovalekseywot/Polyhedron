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
    // На стадии инициализации расставляем индексы ячеек, и начальные потенциалы
    // некоторым образом (u1).
    for (auto cell: *mesh->cells())
    {
        JacobiCellData data;
        data.vol = figure.is_inside(cell->center());
        data.u1 = -V0 * cell->center()[0];
        data.u2 = -V0 * cell->center()[0];
        data.p = 0.0;
        // data.idx = ??
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

void Jacobi::JacobiStage(const NodeList::Part& cells) const {
    // Во время расчетного шага считаем потенциалы на следующей итерации u2,
    // затем записываем u2 на место u1, считаем градиент u1 методом МНК, то
    // есть получаем компоненты Vx, Vy, Vz
    // mu = S_ab/abs(center_b-center_a)
    // u_a * sum(mu) - sum(mu * u_b) = граничное условие или 0
    for (auto cell: cells)
    {
        if (get_state(cell).vol > 0.0)
        {
            auto data = get_state(cell);
            data.u1 = 0.0;
            data.u2 = 0.0;
            data.p = 0.0;
            data.v = Vector3d::Zero();

            int count = 0;
            for (auto &face: cell->faces())
            {
                auto neib = face->neighbor(cell);
                auto neib_data = get_state(neib);

                if (neib_data.vol > 0.0)
                    continue;

                // Соседняя ячейка в жидкости
                count++;
                data.u1 += neib_data.u1;
                data.u2 += neib_data.u2;
                data.v += neib_data.v;
                data.p += 0.5 * neib_data.v.norm() - neib_data.v[0] * V0;
            }
            if (count > 0)
            {
                data.u1 /= count;
                data.u2 /= count;
                data.p /= count;
                data.v /= count;
            }

            set_state(cell, data);
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
                auto data = get_state(neib);
                double mu = face->area() / (neib->center() - cell->center()).dot(face->normal(neib));
                mu_sum += mu;
                neib_sum += data.u1 * mu;
            }
        }

        JacobiCellData data = get_state(cell);
        data.u2 = (cond_sum + neib_sum) / mu_sum;

        set_state(cell, data);
    }
}

void Jacobi::VelocityStage(const NodeList::Part &cells) const {
    // расчёт скорости
    for (auto cell: cells)
    {
        if (get_state(cell).vol > 1e-5 * figure.getMLength())
            continue;
        Matrix3d a = Matrix3d::Zero();
        Vector3d f = Vector3d::Zero();
        JacobiCellData self_data = get_state(cell);

        for (auto &face: cell->faces())
        {
            if (face->flag() != FaceFlag::ORDER || face->flag() == FaceFlag::WALL)
                continue;

            auto data = get_state(face->neighbor(cell));
            Vector3d c = face->neighbor(cell)->center() - cell->center();
            double w = 1.0 / c.squaredNorm();

            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    a(i, j) += w * c(i) * c(j);

            f += w * (self_data.u2 - data.u2) * c;
        }

        _2D(a(2, 2) = 1.0;)
        // TODO     Size<=1) || (Size>4) || (extract_data(src.nestedExpression())!=extract_data(dst)))
        //              && "Aliasing problem detected in inverse(), you need to do inverse().eval() here."
        a = a.inverse().eval();

        self_data.v = a * f;
        set_state(cell, self_data);
    }
}

std::array<double, 2> Jacobi::ErrorsStage(const NodeList::Part &cells) const {
    // расчёт погрешностей
    double eps = 0.0;
    double delta = 0.0;

    for (auto cell: cells)
    {
        if (get_state(cell).vol > 1e-5)
            continue;
        double neib_sum = 0.0, mu_sum = 0.0;
        double cond_sum = 0.0;
        for (auto &face: cell->faces())
        {
            if (face->flag() != FaceFlag::ORDER)
            {
                cond_sum += face->area() * boundary_function(face->flag(), face->center(), -face->normal(cell));
            } else
            {
                auto data = get_state(face->neighbor(cell));
                double mu = face->area() / (face->neighbor(cell)->center() - cell->center()).norm();
                mu_sum += mu;
                neib_sum += data.u2 * mu;
            }
        }

        auto data = get_state(cell);
        eps = max(eps, abs(data.u2 - data.u1));
        delta = max(delta, abs(cond_sum + neib_sum - mu_sum * data.u2));
    }

    return {eps, delta};
}

void Jacobi::UpdateStage(const NodeList::Part &cells) const {
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
    if (first_step) // расстановка cтенок фигуры при первом шаге
    {
        for (auto cell: *mesh->cells())
        {
            JacobiCellData self_data = get_state(cell);
            for (auto &face: cell->faces())
            {
                if (face->flag() != FaceFlag::ORDER)
                {
                    continue;
                }

                auto neib_data = get_state(face->neighbor(cell));
                if (self_data.vol != neib_data.vol)
                {
                    face->set_flag(FaceFlag::WALL);
                }
            }
        }
        first_step = false;
    }

    using std::placeholders::_1;

    mesh->map(std::bind(&Jacobi::JacobiStage, this, _1));

    mesh->map(std::bind(&Jacobi::VelocityStage, this, _1));

    std::vector<double> eps(mesh->n_chunks());
    std::vector<double> delta(mesh->n_chunks());
    mesh->map([this, &eps, &delta](const NodeList::Part& cells) {
       auto res = ErrorsStage(cells);
       eps[cells.part_id()]   = res[0];
       delta[cells.part_id()] = res[1];
    });
    m_eps   = *std::max_element(eps.begin(), eps.end());
    m_delta = *std::max_element(delta.begin(), delta.end());

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
        set_state(child, pdata);
    }
}

void Jacobi::coarse_data(const shared_ptr<Cell> &parent, const vector<shared_ptr<Cell>> &children)
{
    auto chdata = get_state(children[0]);
    set_state(parent, chdata);
}

void Jacobi::print_info(const char *tab) const
{
    std::cout << tab << std::scientific << std::setprecision(4)
              << "Eps = " << m_eps << ", "
              << "Delta = " << m_delta << "\n";
}

double Jacobi::get_cell_param(const shared_ptr<Cell> &cell, const string &name) const
{
    JacobiCellData data = get_state(cell);
    if (name == "u_old")
        return data.u1;
    else if (name == "u_new")
        return data.u2;
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
    else
        throw std::runtime_error("Unknown parameter '" + name + "'");
}

double Jacobi::get_integral_param(const string &name) const
{
    if ("eps") {
        return m_eps;
    } else if ("delta") {
        return m_delta;
    } else {
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
