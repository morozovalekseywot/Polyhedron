#include <utils/surface/surface.hpp>

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE surface
#define BOOST_TEST_NO_LIB

#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(Surface_is_inside)

    BOOST_AUTO_TEST_CASE(is_inside_tetra)
    {
        surf::Surface surf("../examples/figure/tetra.stl");

        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.359, -0.933, -0.479}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.56, -1.30, -0.39}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.00, -1.91, -0.466}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.41, -1.22, -0.21}), true);
    }

    BOOST_AUTO_TEST_CASE(is_inside_prizm)
    {
        surf::Surface surf("../examples/figure/prizm.stl");

        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.55, 0.63, 0.55}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.58, 0.36, 0.34}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.33, 0.15, 0.19}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.99, 0.29, 0.59}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.51, 0.08, 0.04}), true);
    }

    BOOST_AUTO_TEST_CASE(is_inside_pyramid)
    {
        surf::Surface surf("../examples/figure/tetra.stl");

        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.37, 0.84, 0.02}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.67, 0.41, 0.08}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.05, 0.95, 0.86}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.57, 1.57, 0.199}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.27, 0.96, 0.126}), true);
    }

    BOOST_AUTO_TEST_CASE(is_inside_hexahedra)
    {
        surf::Surface surf("../examples/figure/hexahedra.stl");

        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.5, 0, 0}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.29, -0.06, -0.21}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.116, -0.38, -0.45}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.38, 0.43, -0.396}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.165, 0.14, 0.44}), true);
    }

BOOST_AUTO_TEST_SUITE_END()