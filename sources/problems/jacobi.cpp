#include <control/configuration.h>
#include <core/face/face.h>
#include <core/cell/cell.h>
#include <core/cell/data_holder.h>
#include <core/mesh/mesh.h>
#include <core/generator/geometry_generator.h>
#include <core/generator/rectangle_generator.h>
#include <problems/jacobi.h>

Jacobi::Jacobi(const Configuration &config)
{
    figure = surf::Surface("examples/figure/pyramid2.stl");
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
        data.u1 = 0.0;
        data.u2 = 0.0;
        // data.idx = ??
        set_state(cell, data);
    }
}

double boundary_function(FaceFlag flag,const Vector3d &vec, const Vector3d &n)
{
    if(flag == FaceFlag::INFLOW)
        return 10.0;

    if(flag == FaceFlag::OUTFLOW)
        return -10.0;

    return 0.0;
}

double Jacobi::solution_step(Mesh *mesh)
{
    if(first_step) // расстановка cтенок фигуры при первом шаге
    {
        for (auto cell: *mesh->cells()) {
            JacobiCellData self_data = get_state(cell);
            for (auto &face: cell->faces()) {
                if (face->flag() != FaceFlag::ORDER) {
                    continue;
                }

                auto neib_data = get_state(face->neighbor(cell));
                if (self_data.vol != neib_data.vol) {
                    face->set_flag(FaceFlag::WALL);
                }
            }
        }
        first_step = false;
    }

    // Во время расчетного шага считаем потенциалы на следующей итерации u2,
    // затем записываем u2 на место u1, считаем градиент u1 методом МНК, то
    // есть получаем компоненты Vx, Vy, Vz
    // mu = S_ab/abs(center_b-center_a)
    // u_a * sum(mu) - sum(mu * u_b) = граничное условие или 0
    for (auto cell: *mesh->cells())
    {
        if (get_state(cell).vol > 0)
        {
            auto data = get_state(cell);
            data.u2 = 0.0;
            data.v = Vector3d::Zero();
            set_state(cell, data);
            continue;
        }

        double neib_sum = 0.0, mu_sum = 0.0;
        double cond_sum = 0.0;
        for (auto &face: cell->faces())
        {
//            if(face->flag() == FaceFlag::WALL)
//            {
//                mu_sum += face->area() / (face->neighbor(cell)->center() - cell->center()).norm();
//                continue;
//            }

            if (face->flag() != FaceFlag::ORDER)
            {
                cond_sum += face->area() * boundary_function(face->flag(),face->center(), -face->normal(cell));
            } else
            {
                auto data = get_state(face->neighbor(cell));
                double mu = face->area() / (face->neighbor(cell)->center() - cell->center()).norm();
                mu_sum += mu;
                neib_sum += data.u1 * mu;
            }
        }

        JacobiCellData data = get_state(cell);
        data.u2 = (cond_sum + neib_sum) / mu_sum;

        set_state(cell, data);
    }

    // расчёт скорости
    for (auto cell: *mesh->cells())
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

    // расчёт погрешностей
    m_eps = 0.0;
    m_delta = 0.0;
    for (auto cell: *mesh->cells())
    {
        if (get_state(cell).vol > 1e-5)
            continue;
        double neib_sum = 0.0, mu_sum = 0.0;
        double cond_sum = 0.0;
        for (auto &face: cell->faces())
        {
            if (face->flag() != FaceFlag::ORDER)
            {
                cond_sum += face->area() * boundary_function(face->flag(),face->center(), -face->normal(cell));
            } else
            {
                auto data = get_state(face->neighbor(cell));
                double mu = face->area() / (face->neighbor(cell)->center() - cell->center()).norm();
                mu_sum += mu;
                neib_sum += data.u2 * mu;
            }
        }

        auto data = get_state(cell);
        m_eps = max(m_eps, abs(data.u2 - data.u1));
        m_delta = max(m_delta, abs(cond_sum + neib_sum - mu_sum * data.u2));
    }

    for (auto cell: *mesh->cells())
    {
        auto data = get_state(cell);
        data.u1 = data.u2;
        data.u2 = 0.0;
        set_state(cell, data);
    }

    m_time += 1.0;
    return 0;
}

AdaptationFlag Jacobi::adaptation_criterion(const shared_ptr<Cell> &cell)
{
    if (m_time > 0.0) {
        return AdaptationFlag::NONE;
    }

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
    std::cout << std::setprecision(4) << "Eps = " << m_eps;
    std::cout << ", Delta = " << m_delta << "\n";
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
    else
        throw std::runtime_error("Unknown parameter '" + name + "'");
}

double Jacobi::get_integral_param(const string &name) const
{
    return 0;
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
