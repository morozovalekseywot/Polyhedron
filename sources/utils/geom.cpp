//
// Created by 159-mrv on 8/14/18.
//

#include <core/vertex/vertex.h>
#include <core/cell/cell.h>
#include <utils/geom.h>


Vector2d geom::to_polar(const Vector2d &vertex) {
    return { vertex.norm(), atan2(vertex[1], vertex[0]) };
}

Vector2d geom::to_cartesian(double r, double phi) {
    return { r * cos(phi), r * sin(phi) };
}

Vector2d geom::lines_intersection(const Vector2d &v1, const Vector2d &v2,
                                  const Vector2d &v3, const Vector2d &v4) {
    Vector2d a = v2 - v1;
    Vector2d b = v4 - v3;

    double D = a[0] * b[1] - a[1] * b[0];
    double t = ((v3[0] - v1[0]) * b[1] - (v3[1] - v1[1]) * b[0]) / D;

    return v1 + t * a;
}

double geom::sqr_distance(const Vector3d &v, const Vector3d &v1, const Vector3d &v2) {
    Vector3d a = v2 - v1;
    double xi = (v - v1).dot(a) / a.squaredNorm();

    if (xi > 1.0) {
        xi = 1.0;
    }
    else if (xi < 0.0) {
        xi = 0.0;
    }

    Vector3d v_int = v1 + a * xi;
    return (v - v_int).squaredNorm();
}

double geom::sqr_distance(const Vector3d &v, const Vector3d &v1, const Vector3d &v2,
                          const Vector3d &v3, const Vector3d &v4) {
/*
    Vector3d c = 0.25 * (v1 + v2 + v3 + v4);
    Vector3d n = geom::normal({v1, v2, v3, v4});

    double h = (v - c).dot(n);
    Vector3d vt = v - h * n;

    bool inside = (v1 - c).cross(v2 - c).dot(n) > 0.0 &&
                  (v2 - c).cross(v3 - c).dot(n) > 0.0 &&
                  (v3 - c).cross(v4 - c).dot(n) > 0.0 &&
                  (v4 - c).cross(v1 - c).dot(n) > 0.0;

    if (inside) {
        return h * h;
    }
    */

    double d1 = geom::sqr_distance(v, v1, v2);
    double d2 = geom::sqr_distance(v, v2, v3);
    double d3 = geom::sqr_distance(v, v3, v4);
    double d4 = geom::sqr_distance(v, v4, v1);
    double d5 = geom::sqr_distance(v, v1, v3);
    double d6 = geom::sqr_distance(v, v2, v4);
    return min(min(d1, d2), min(min(d3, d4), min(d5, d6)));
}

double geom::polygon_area(const vector<Vector2d>& vertices) {
    double area = 0.0;
    size_t j = vertices.size() - 1;

    for (size_t i = 0; i < vertices.size(); ++i) {
        area += (vertices[j][0] + vertices[i][0]) *
                (vertices[j][1] - vertices[i][1]);
        j = i;
    }
    return -area / 2;
}

double geom::polygon_area(const array<Vector2d, 3>& v) {
    return 0.5 * fabs(
            (v[0][0] - v[2][0]) * (v[1][1] - v[2][1]) -
            (v[1][0] - v[2][0]) * (v[0][1] - v[2][1])
    );
}

double geom::polygon_area(const array<Vector3d, 3>& v) {
    Vector3d a = v[2] - v[0];
    Vector3d b = v[1] - v[0];
    return 0.5 * a.cross(b).norm();
}

double geom::polygon_area(const array<Vector2d, 4>& v) {
    return 0.5 * fabs(
            (v[0][0] - v[1][0]) * (v[0][1] + v[1][1]) +
            (v[1][0] - v[2][0]) * (v[1][1] + v[2][1]) +
            (v[2][0] - v[3][0]) * (v[2][1] + v[3][1]) +
            (v[3][0] - v[0][0]) * (v[3][1] + v[0][1])
    );
}

double geom::polygon_area(const array<Vector3d, 4>& v) {
    Vector3d a = v[1] - v[0];
    Vector3d b = v[2] - v[1];
    Vector3d c = v[3] - v[2];
    Vector3d d = v[0] - v[3];
    return 0.25 * (a.cross(b).norm() + b.cross(c).norm() +
                   c.cross(d).norm() + d.cross(a).norm());
}

Vector3d geom::normal(const array<Vector3d, 4>& v) {
    Vector3d a = v[1] - v[0];
    Vector3d b = v[2] - v[1];
    Vector3d c = v[3] - v[2];
    Vector3d d = v[0] - v[3];

    Vector3d norm = a.cross(b) + b.cross(c) + c.cross(d) + d.cross(a);
    norm.normalize();
    return norm;
}

Vector2d geom::symmetric_point(const Vector2d& p, const Vector2d& v1, const Vector2d& v2) {
    // Нормаль к прямой
    Vector2d normal = (v2 - v1).normalized();
    normal = { normal[1], -normal[0] };

    Vector2d rc = v1 - p;

    //rc_n = (rc,n)*n - vector along normal_2d
    double rc_n_proj = rc.dot(normal);
    Vector2d rc_n = normal * rc_n_proj;

    return p + 2.0 * rc_n;
}

double geom::split_cell(shared_ptr<Cell> cell, const Vector2d& v1, const Vector2d& v2) {
    // Вершины ячейки
    Vector2d lbv = cell->left_bottom_vertex()->v_2d();
    Vector2d ltv = cell->left_top_vertex()->v_2d();
    Vector2d rbv = cell->right_bottom_vertex()->v_2d();
    Vector2d rtv = cell->right_top_vertex()->v_2d();

    // Полярные границы ячейки
    double r1 = lbv.norm();
    double r2 = rtv.norm();

    Vector2d inter_top    = lines_intersection(v1, v2, ltv, rtv);
    Vector2d inter_bottom = lines_intersection(v1, v2, lbv, rbv);
    Vector2d inter_left   = lines_intersection(v1, v2, lbv, ltv);
    Vector2d inter_right  = lines_intersection(v1, v2, rbv, rtv);

    double r_left  = inter_left.norm();
    double r_right = inter_right.norm();

    double full_area = cell->volume();
    if (r_left < r1) {
        if (r_right < r1) {
            // Линия проходит ниже ячейки
            return 1.0;
        }
        else if (r_right < r2) {
            // Отсекается правый нижний угол
            array<Vector2d, 3> tri = {inter_bottom, rbv, inter_right};
            double part_area = geom::polygon_area(tri);
            return (full_area - part_area) / full_area;
        }
        else {
            // маловероятный случай, линия прохоит вертикально
            array<Vector2d, 4> quad = {inter_bottom, rbv, rtv, inter_top};
            double part_area = geom::polygon_area(quad);
            return (full_area - part_area) / full_area;
        }
    }
    else if (r_left < r2) {
        if (r_right < r1) {
            // Отсекается левый нижний угол
            array<Vector2d, 3> tri = {lbv, inter_bottom, inter_left};
            double part_area = polygon_area(tri);
            return (full_area - part_area) / full_area;
        }
        else if (r_right < r2) {
            // Отсекается нижняя часть
            array<Vector2d, 4> quad = {lbv, rbv, inter_right, inter_left};
            double part_area = geom::polygon_area(quad);
            return (full_area - part_area) / full_area;
        }
        else {
            // Отсекается левый верхний угол
            array<Vector2d, 3> tri = {inter_left, inter_top, ltv};
            double part_area = geom::polygon_area(tri);
            return part_area / full_area;
        }
    }
    else {
        if (r_right < r1) {
            // Маловероятный случай, линия проходит вертикально
            array<Vector2d, 4> quad = {inter_bottom, rbv, rtv, inter_top};
            double part_area = geom::polygon_area(quad);
            return part_area / full_area;
        } else if (r_right < r2) {
            // Отсекается правый верхний угол
            array<Vector2d, 3> tri = {inter_right, rtv, inter_top};
            double part_area = geom::polygon_area(tri);
            return part_area / full_area;
        } else {
            // Линия проходит над ячейкой
            return 0.0;
        }
    }
};