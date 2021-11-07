#include <control/configuration.h>
#include <core/face/face.h>
#include <core/cell/cell.h>
#include <core/cell/data_holder.h>
#include <core/mesh/mesh.h>
#include <core/generator/geometry_generator.h>
#include <core/generator/rectangle_generator.h>
#include <problems/jacobi.h>

Jacobi::Jacobi(const Configuration &config)
{}


uint Jacobi::cell_data_size() const
{
    return static_cast<uint>(ceil(double(sizeof(JacobiCellData)) / sizeof(double)));
}

/*
double is_inside(const Vector3d &r)
{
    // задаём трубу
    Vector3d right_up = {0.7, 0.2, 0.0};
    Vector3d right_down = {0.7, -0.2, 0.0};
    Vector3d left_up = {-0.70, 0.2, 0.0};
    Vector3d left_down = {-0.7, 0.2, 0.0};

    Vector3d up = right_up - left_up;
    Vector3d right = right_down - right_up;
    Vector3d down = left_down - right_down;
    Vector3d left = left_up - left_down;

    std::vector<std::pair<Vector3d, Vector3d>> sides(4);// пара: (нормаль, центр)
    sides[0] = std::make_pair(Vector3d{up.x(), up.y(), 0.0}, up / 2);
    sides[1] = std::make_pair(Vector3d{right.x(), right.y(), 0.0}, right / 2);
    sides[2] = std::make_pair(Vector3d{down.x(), down.y(), 0.0}, down / 2);
    sides[3] = std::make_pair(Vector3d{left.x(), left.y(), 0.0}, left / 2);
    for (auto &side: sides)
    {
        if (side.first.x() * side.second.x() + side.first.y() * side.second.y() +
            side.first.z() * side.second.z() < 0.0)
        {
            side.first = -side.first;
        }
    }

    for (auto &side: sides)
    {
        Vector3d vec = side.second - r;
        if (side.first.x() * vec.x() + side.first.y() * vec.y() +
            side.first.z() * vec.z() < 0.0)
            return 0.0;
    }

    return 1.0;
}
*/

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

//        for (auto &face: cell->faces())
//        {
//            auto neib = face->neighbor(cell);
//
//            if (figure.is_inside(neib->center()) == 0.0)
//            {
//                // Сколько ставить отличие по y??
//                Vector3d r = neib->center() - cell->center();
//                if (r.x() < 0.0 && abs(r.y()) < abs(r.x()) / 20)
//                {
//                    data.v = {10.0, 0.0, 0.0};
//                    // data.condition = 10.0; ??
//                } else if (r.x() > 0.0 && abs(r.y()) < abs(r.x()) / 20)
//                {
//                    data.v = {20.0, 0.0, 0.0};
//                    // data.condition = 20.0; ??
//                }
//                break;
//            }
//        }
//        data.u1 = double(rand()) / RAND_MAX;
//        data.u2 = 0.0;
        set_state(cell, data);
    }
}

double boundary_function(const Vector3d &vec, const Vector3d &n)
{
    if (vec.x() < 0 && n.x() < -0.5)
    {
        return 10.0;
    }
    if (vec.x() > 0 && vec.y() < 0 && n.x() > 0.5)
    {
        return -20.0;
    }

    return 0;
}

double Jacobi::solution_step(Mesh *mesh)
{
    // Во время расчетного шага считаем потенциалы на следующей итерации u2,
    // затем записываем u2 на место u1, считаем градиент u1 методом МНК, то
    // есть получаем компоненты Vx, Vy, Vz
    // mu = S_ab/abs(center_b-center_a)
    // u_a * sum(mu) - sum(mu * u_b) = граничное условие или 0
    for (auto cell: *mesh->cells())
    {
        double neib_sum = 0.0, mu_sum = 0.0;
        double cond_sum = 0.0;
        for (auto &face: cell->faces())
        {
            if (face->flag() != FaceFlag::ORDER)
            {
                cond_sum += face->area() * boundary_function(face->center(), -face->normal(cell));
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
        Matrix3d a = Matrix3d::Zero();
        Vector3d f = Vector3d::Zero();
        JacobiCellData self_data = get_state(cell);

        for (auto &face: cell->faces())
        {
            if (face->flag() != FaceFlag::ORDER)
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
        double neib_sum = 0.0, mu_sum = 0.0;
        double cond_sum = 0.0;
        for (auto &face: cell->faces())
        {
            if (face->flag() != FaceFlag::ORDER)
            {
                cond_sum += face->area() * boundary_function(face->center(), -face->normal(cell));
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
    std::cout << ",Delta = " << m_delta << "\n";
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
