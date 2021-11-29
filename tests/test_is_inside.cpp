#include <utils/surface/surface.hpp>

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE surface
#define BOOST_TEST_NO_LIB

#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(Surface_is_inside)

    BOOST_AUTO_TEST_CASE(is_inside_tetra)
    {
        surf::Surface surf("../examples/figure/tetra.stl");

        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.359, -0.933, -0.479}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.56, -1.30, -0.39}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.00, -1.91, -0.466}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.41, -1.22, -0.21}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.5, 0, 0}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.5, 0, 0}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1.1, -1.54, 0.36}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-1, 1.7, 0.82}), false);
    }

    BOOST_AUTO_TEST_CASE(is_inside_prizm)
    {
        surf::Surface surf("../examples/figure/prizm.stl");

        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.55, 0.63, 0.55}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.58, 0.36, 0.34}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.33, 0.15, 0.19}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.99, 0.29, 0.59}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.51, 0.08, 0.04}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.5, 0, 0}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.16, 0.99, 0.77}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{2.48, 0.48, 0.43}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.94, 0.62, -0.36}), false);
    }

    BOOST_AUTO_TEST_CASE(is_inside_pyramid)
    {
        surf::Surface surf("../examples/figure/pyramid.stl");

        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.37, 0.84, 0.02}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.67, 0.41, 0.08}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.05, 0.95, 0.86}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.57, 1.57, 0.199}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.27, 0.96, 0.126}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.38, 1.81, 0.37}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.37, -0.53, 0.4}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.8, 1.98, 0.7}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.5, -0.25, 0.087}), false);
    }

    BOOST_AUTO_TEST_CASE(is_inside_hexahedra)
    {
        surf::Surface surf("../examples/figure/hexahedra.stl");

        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.5, 0, 0}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.29, -0.06, -0.21}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.116, -0.38, -0.45}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.38, 0.43, -0.396}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.165, 0.14, 0.44}), true);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.2, 0.147, -0.798}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{-0.3, 0.25, 0.2}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{0.48, -0.08, -0.89}), false);
        BOOST_CHECK_EQUAL(surf.is_inside(Vector3d{1.36, -0.16, -0.26}), false);
    }

BOOST_AUTO_TEST_SUITE_END()