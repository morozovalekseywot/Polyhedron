#include <utils/surface/surface.hpp>


surf::Surface::Surface(const std::string &filename)
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
//            std::vector<std::array<Vector3d, 4>> raw_triangles;
    std::string t;
    getline(file, t); // игнорируем имя solid
    while (!file.eof())
    {
        file >> t; // facet
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

    /// Вариант 1. Тупой, за O(N^2).
    /// Проходим по треугольникам в raw_triangles, треугольник кидаем
    /// в m_triangles, затем у каждого треугольника проходим по вершинам.
    /// Сначала ищем вершину в m_vertices (простым перебором), если
    /// вершины там нет (с точностью до eps*m_length), тогда добавляем
    /// новую вершину в конец m_vertices, а индекс этой вершины
    /// записываем в вершину к m_triangles.back().

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
        m_triangles[size] = Triangle(size, {idx1, idx2, idx3}, face.normal);
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

double surf::Surface::getMLength() const
{
    return m_length;
}

const std::vector<surf::Vertex> &surf::Surface::getMVertices() const
{
    return m_vertices;
}

const std::vector<surf::Triangle> &surf::Surface::getMTriangles() const
{
    return m_triangles;
}
