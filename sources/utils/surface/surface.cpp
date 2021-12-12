#include <utils/surface/surface.hpp>

namespace surf {

Surface::Surface(const std::string &filename)
{
    /// В конструкторе необходимо заполнить массивы m_vertices
    /// и m_triangles.

    /// Считываем треугольники из файла в промежуточный массив,
    /// где каждый треугольник представляется набором из трех
    /// вершин v1, v2, v3 и нормали n, далее работаем с ним
    std::string path = "../" + filename;
    std::ifstream file(path, std::ios_base::in);
    if (!file.is_open()) // если файл не был открыт
    {
        throw std::runtime_error("file isn't open\n");
    }
    std::vector<Face> faces;
    std::string t;
    getline(file, t); // игнорируем имя solid
    while (!file.eof())
    {
        file >> t; // facet
        if (t == "endsolid")
            break;
        file >> t; // normal
        Vector3d n;
        file >> n.x() >> n.y() >> n.z();
        getline(file, t); // ignore /n
        getline(file, t); // outer loop
        std::array<Vector3d, 3> verts;
        for (auto &v: verts)
        {
            file >> t; // vertex
            file >> v.x() >> v.y() >> v.z();
        }
        faces.emplace_back(Face(verts, n));
        getline(file, t); // ignore /n
        getline(file, t); // endloop
        getline(file, t); // endfacet
    }
    file.close(); // закрываем файл

    /// По ходу считывания из файла вычисляем m_length. Более точно:
    /// строим BoundingBox, затем m_length считаем как длину диагонали.
    /// В дальнейшем m_length это характерный размер.
    BBox bbox(faces);
    m_length = bbox.len.norm();
    /// Далее самая сложная операция, удалить дубликаты вершин.

    /// Вариант 2. Адекватный, за O(N log(N)).
    /// Вершины добавляются параллельно в m_vertices и в некоторый
    /// set, для которого пишется определенный компаратор.
    /// Вместо прохода по всем добавленным вершинам, новая вершина
    /// просто ищется в set.
    /// Компаратор следующий: вершины эквивалентны, если совпадают
    /// с точностью до eps*m_length, иначе сравниваются их проекции
    /// на некоторую довольно случайную ось.
    /// Хотя у меня есть сомнения, что вы осознаете, что тут написано.
    /// Пусть eps = 1e-8

    auto comp = [this](const Eigen::Vector3d &a, const Eigen::Vector3d &b) -> bool
    {
        if (equal_d(a.x(), b.x()))
        {
            if (equal_d(a.y(), b.y()))
            {
                if (equal_d(a.z(), b.z()))
                    return false;
                else
                    return a.z() < b.z();
            } else
                return a.y() < b.y();
        } else
            return a.x() < b.x();
    };

    std::set<Eigen::Vector3d, decltype(comp)> set_vtx(comp);
    for (auto &face: faces)
    {
        set_vtx.insert(face.vects[0]);
        set_vtx.insert(face.vects[1]);
        set_vtx.insert(face.vects[2]);
    }

    m_vertices.resize(set_vtx.size());
    size_t size = 0;
    for (auto &v: set_vtx)
        m_vertices[size++] = Vertex(v);

    size = 0;
    m_triangles.resize(faces.size());
    for (auto &face: faces)
    {
        // ищем индексы наших вершин в массиве m_vertices бин поиском (lower_bound)
        auto it1 = std::lower_bound(m_vertices.begin(), m_vertices.end(), face.vects[0], comp);
        auto it2 = std::lower_bound(m_vertices.begin(), m_vertices.end(), face.vects[1], comp);
        auto it3 = std::lower_bound(m_vertices.begin(), m_vertices.end(), face.vects[2], comp);
        size_t idx1 = std::distance(m_vertices.begin(), it1);
        size_t idx2 = std::distance(m_vertices.begin(), it2);
        size_t idx3 = std::distance(m_vertices.begin(), it3);
        m_triangles[size] = Triangle(size, {idx1, idx2, idx3}, face.normal, face.area);
        size++;
    }
    /// Ладно, пофиг, реализовали первый варинат. Теперь у нас есть
    /// m_vertices только с положениями, без дубликатов и
    /// m_triangles, которые на них ссылаются.

    /// Проходим по треугольникам в m_triangles, у них по вершинам,
    /// для вершин записываем смежные треугольники.
    for (int i = 0; i < m_triangles.size(); ++i)
    {
        std::array<size_t, 3> verts = m_triangles[i].vertices;
        for (int j = 0; j < 3; j++)
            m_vertices[verts[j]].set_next_triangle(i);
    }

    /// Проходим по треугольникам в m_triangles, проходим по их
    /// вершинам, у вершин в треугольниках ищем смежные треугольники.
    for (auto &triangle: m_triangles)
    {
        size = 0;
        // создаём вектор с индексами треугольников, куда будем сохранять все индексы треугольников у вершин
        std::vector<size_t> list_of_triangles(3 * m_vertices[triangle.vertices[0]].max_triangles);
        // копируем все индексы треугольников из вершин треугольника
        for (int idx = 0; idx < 3; idx++)
        {
            size_t verts_idx = triangle.vertices[idx];
            std::copy(m_vertices[verts_idx].triangles.begin(), m_vertices[verts_idx].triangles.end(), list_of_triangles.begin() + size);
            size += m_vertices[verts_idx].triangles.size();
        }

        list_of_triangles.resize(size);
        // сортируем и ищем индексы которые встречаются два раза, и надо исключить номер самого треугольника (встречается 3 раза)
        std::sort(list_of_triangles.begin(), list_of_triangles.end());
        size = 0;
        for (int i = 0; i < list_of_triangles.size() - 1; i++)
        {
            if (list_of_triangles[i] == list_of_triangles[i + 1] && list_of_triangles[i] != triangle.index)
                triangle.triangles[size++] = list_of_triangles[i];
        }
    }
}

double Surface::getMLength() const
{
    return m_length;
}

const std::vector<Vertex> &Surface::getMVertices() const
{
    return m_vertices;
}

const std::vector<Triangle> &Surface::getMTriangles() const
{
    return m_triangles;
}

bool Surface::is_inside(const Vector3d &v_) const
{
    /// Здесь простая реализация. Будем использовать простейшее правило
    /// четный-нечетный (even-odd rule). Из какой-нибудь точки выпускается
    /// луч в произвольном направлении. Затем отыскиваются пересечения луча
    /// и треугольников поверхности.
    Vector3d v = v_;
    Vector3d s = {0.257381649640262, 0.5428934530125299, 0.7993881317011088};
    s *= 1000 * m_length;
//    for (auto &idx: m_triangles[m_vertices.back().triangles[0]].vertices) // берём треугольник смежный последней точке
//        s = s + m_vertices[idx].v;
//    s /= 3.0;
//    s = v + 100 * (s - v); // продлил луч
    Vector3d l = s - v; // направляющий вектор прямой

    int count = 0; // количество пересечений

    auto scalar = [](const Vector3d &a, const Vector3d &b) -> double
    {
        return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
    };

    for (auto &triangle: m_triangles)
    {
        // В первую очередь, проверяется, пересекает ли контрольный луч плоскость просматриваемого
        // треугольника: точки v и 𝑆 должны лежать с разных сторон или скалярные произведения
        // нормали с радиус-векторами от любой вершины до этих точек должны иметь разные знаки. В
        // противном случае – переход к следующему треугольнику без дальнейших проверок.
        Eigen::Vector3d first = m_vertices[triangle.vertices[0]].v;
        if (scalar(triangle.normal, v - first) * scalar(triangle.normal, s - first) > 0.0)
            continue;

        if (abs(scalar(triangle.normal, first - v)) / (first - v).norm() < 1e-5) // точка лежит в плоскости треугольника
        {
            std::array<Vector3d, 3> arr;
            for (int i = 0; i < 3; i++)
                arr[i] = m_vertices[triangle.vertices[i]].v;

            double area_ = ((v - arr[0]).cross(v - arr[1])).norm() + ((v - arr[0]).cross(v - arr[2])).norm() + ((v - arr[1]).cross(v - arr[2])).norm();
            return (abs(area_ - 2 * triangle.area) < triangle.area * 1e-5);
        }
        // ур-е плоскости: Ax + By + Cz + D = 0
        // ур-е прямой: x = v.x + m * t, y = v.y + p * t, z = v.z + l * t
        // A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2)
        // B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2)
        // C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2)
        // -D = x1 (y2 z3 - y3 z2) + x2 (y3 z1 - y1 z3) + x3 (y1 z2 - y2 z1)
        std::array<Vector3d, 3> arr;
        for (int i = 0; i < 3; i++)
            arr[i] = m_vertices[triangle.vertices[i]].v;

        double A, B, C, D;
        A = arr[0].y() * (arr[1].z() - arr[2].z()) + arr[1].y() * (arr[2].z() - arr[0].z()) + arr[2].y() * (arr[0].z() - arr[1].z());
        B = arr[0].z() * (arr[1].x() - arr[2].x()) + arr[1].z() * (arr[2].x() - arr[0].x()) + arr[2].z() * (arr[0].x() - arr[1].x());
        C = arr[0].x() * (arr[1].y() - arr[2].y()) + arr[1].x() * (arr[2].y() - arr[0].y()) + arr[2].x() * (arr[0].y() - arr[1].y());
//                D = -(arr[0].x() * (arr[1].y() * arr[2].z() - arr[1].z() * arr[2].y()) + arr[1].x() * (arr[2].y() * arr[0].z() - arr[2].z() * arr[0].y()) +
//                      arr[2].x() * (arr[0].y() * arr[1].z() - arr[0].z() * arr[1].y()));
        D = -A * arr[0].x() - B * arr[0].y() - C * arr[0].z();

        double t = -(D + A * v.x() + B * v.y() + C * v.z()) / (A * l.x() + B * l.y() + C * l.z());
        if (t < 0.0)
            continue;
        Vector3d x = {v.x() + t * l.x(), v.y() + t * l.y(), v.z() + t * l.z()}; // точка пересечения с плоскостью треугольника
        std::array<double, 3> arr_area = {((x - arr[0]).cross(x - arr[1])).norm(), ((x - arr[0]).cross(x - arr[2])).norm(),
                                          ((x - arr[1]).cross(x - arr[2])).norm()};
        //double area_ = ((x - arr[0]).cross(x - arr[1])).norm() + ((x - arr[0]).cross(x - arr[2])).norm() + ((x - arr[1]).cross(x - arr[2])).norm();
        if (abs(arr_area[0] + arr_area[1] + arr_area[2] - 2 * triangle.area) > triangle.area * 1e-5)
            continue; // точка не в треугольнике
        else
        {
            for (auto &area: arr_area)
            {
                if (area < triangle.area * 1e-5)
                {
                    v = (250 * v + (arr[0] + arr[1] + arr[2]) / 3) / 251; // сдвигаем точку внутрь треугольника
                    break;
                }
            }
            count++;
        }
    }

    /// Если пересечений нечетное количество,
    /// значит точка лежит внутри поверхности, иначе -- снаружи.
    /// Если точка и луч достаточно случайные, то можно считать, что
    /// вероятность попадания луча на границу треугольника равна нулю
    /// (все должно работать короче).

    return count % 2 != 0;
}

void Surface::centering() {
    Vector3d vmin = m_vertices[0];
    Vector3d vmax = m_vertices[0];
    for (auto &vert: m_vertices) {
        for (int i = 0; i < 3; ++i) {
            vmin[i] = std::min(vmin[i], vert.v[i]);
            vmax[i] = std::max(vmax[i], vert.v[i]);
        }
    }
    Vector3d c = 0.5 * (vmin + vmax);
    for (auto& vert: m_vertices) {
        vert.v -= c;
    }
}

void Surface::rotate(Vector3d n, double phi) {
    // Ось вращения
    n.normalize();

    double nx = n[0];
    double ny = n[1];
    double nz = n[2];

    // Синус и косинус угла вращения
    double c = std::cos(phi);
    double s = std::sin(phi);

    // матрица поворота
    Matrix3d R;
    R(0, 0) = c + (1 - c) * nx * nx;
    R(0, 1) = (1 - c) * nx * ny - s * nz;
    R(0, 2) = (1 - c) * nx * nz + s * ny;

    R(1, 0) = (1 - c) * nx * ny + s * nz;
    R(1, 1) = c + (1 - c) * ny * ny;
    R(1, 2) = (1 - c) * ny * nz - s * nx;

    R(2, 0) = (1 - c) * nx * nz - s * ny;
    R(2, 1) = (1 - c) * ny * nz + s * nx;
    R(2, 2) = c + (1 - c) * nz * nz;

    for (auto &vert: m_vertices) {
        vert.v = R * vert.v;
    }

    for (auto &tri: m_triangles) {
        tri.normal = R * tri.normal;
    }
}

}