#include <utils/surface/surface.hpp>

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE surface
#define BOOST_TEST_NO_LIB

#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using Eigen::Vector3d;
BOOST_AUTO_TEST_SUITE(Surface)

    BOOST_AUTO_TEST_CASE(surface_constructor)
    {
        surf::Surface surf("../examples/pyramid.stl");
        BOOST_CHECK_CLOSE_FRACTION(surf.getMLength(), std::sqrt(2 * 2 + 2 * 2 + 0), 0.0001); // (0,0,0), (2,2,0)
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
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.37 , 0.84 , 0.02}),true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.67 , 0.41 , 0.08 }),true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.05,0.95,0.86}),true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.57,1.57,0.199}),true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.27,0.96,0.126}),true);
    }

BOOST_AUTO_TEST_SUITE_END()