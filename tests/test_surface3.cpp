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
        surf::Surface surf("../examples/tetra.stl");
        BOOST_CHECK_CLOSE_FRACTION(surf.getMLength(), std::sqrt(1.5 * 1.5 + 1 + 1), 0.0001); // (-2,-2,-0.5), (-0.5,-1,0.5)
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
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.359 , -0.933 , -0.479}),true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.56 , -1.30 , -0.39 }),true);

        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.00,-1.91,-0.466}),true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.41,-1.22,-0.21}),true);
    }

BOOST_AUTO_TEST_SUITE_END()