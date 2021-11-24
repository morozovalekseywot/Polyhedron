#include <utils/surface/surface.hpp>

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE surface
#define BOOST_TEST_NO_LIB

#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(Surface)

    BOOST_AUTO_TEST_CASE(surface_constructor)
    {
        surf::Surface surf("../examples/Prizm.stl");
        BOOST_CHECK_CLOSE_FRACTION(surf.getMLength(), std::sqrt(2.5 * 2.5 + 1 + 1), 0.0001); // (0,0,0), (2.5,1,1)
        std::vector<surf::Triangle> triangles = surf.getMTriangles();
        std::vector<surf::Vertex> verts = surf.getMVertices();
        for (auto &triangle: triangles)
        {
            std::array<size_t, 3> t_verts = triangle.vertices;
            for (int i = 0; i < 3; i++)
            {
                auto it = std::find(verts[t_verts[i]].triangles.begin(), verts[t_verts[i]].triangles.end(), triangle.index);
                BOOST_CHECK_EQUAL(*it, triangle.index); // проверяем что в каждой вершине треугольника хранится индекс самого треугольника
            }
        }

        for (auto &triangle: triangles)
        {
            std::array<size_t, 3> trs = triangle.triangles;
            for (int i = 0; i < 3; i++)
            {
                auto it = std::find(triangles[trs[i]].triangles.begin(), triangles[trs[i]].triangles.end(), triangle.index);
                BOOST_CHECK_EQUAL(*it, triangle.index); // проверяем что у всех соседних треугольников наш треугольник записан как сосед
            }
        }
    }

BOOST_AUTO_TEST_SUITE_END()
