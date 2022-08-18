#include "mars_test.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include "mars_bisection.hpp"
#include "mars_lagrange_element.hpp"
#include "mars_mesh.hpp"
#include "mars_quality.hpp"
#include "mars_simplex.hpp"
#include "mars_vtk_writer.hpp"

namespace mars {

    void test_midpoint_index() {
        using namespace mars;

        // tri
        {
            const auto i01 = midpoint_index<2>(0, 1);
            assert(i01 == 3);
            const auto i02 = midpoint_index<2>(0, 2);
            assert(i02 == 4);
            const auto i03 = midpoint_index<2>(1, 2);
            assert(i03 == 5);
        }

        // tet
        {
            const auto i01 = midpoint_index<3>(0, 1);
            assert(i01 == 4);
            const auto i02 = midpoint_index<3>(0, 2);
            assert(i02 == 5);
            const auto i03 = midpoint_index<3>(0, 3);
            assert(i03 == 6);

            const auto i04 = midpoint_index<3>(1, 2);
            assert(i04 == 7);
            const auto i13 = midpoint_index<3>(1, 3);
            assert(i13 == 8);

            const auto i23 = midpoint_index<3>(2, 3);
            assert(i23 == 9);
        }

        // pentatope
        {
            const auto i01 = midpoint_index<4>(0, 1);
            assert(i01 == 5);
            const auto i02 = midpoint_index<4>(0, 2);
            assert(i02 == 6);
            const auto i03 = midpoint_index<4>(0, 3);
            assert(i03 == 7);
            const auto i04 = midpoint_index<4>(0, 4);
            assert(i04 == 8);

            const auto i12 = midpoint_index<4>(1, 2);
            assert(i12 == 9);
            const auto i13 = midpoint_index<4>(1, 3);
            assert(i13 == 10);
            const auto i14 = midpoint_index<4>(1, 4);
            assert(i14 == 11);

            const auto i23 = midpoint_index<4>(2, 3);
            assert(i23 == 12);
            const auto i24 = midpoint_index<4>(2, 4);
            assert(i24 == 13);

            const auto i34 = midpoint_index<4>(3, 4);
            assert(i34 == 14);
        }
    }

    void test_red_refinement_interpolator() {
        using namespace mars;

        SimplexInterpolator<2> mat2;
        red_refinement_interpolator<2>(mat2);
        mat2.describe(std::cout);

        SimplexInterpolator<3> mat3;
        red_refinement_interpolator<3>(mat3);
        mat3.describe(std::cout);

        SimplexInterpolator<4> mat4;
        red_refinement_interpolator<4>(mat4);
        mat4.describe(std::cout);
    }

    void test_red_refinement() {
        using namespace mars;

        Pentatope4 ptope;
        ptope.nodes[0] = 0;
        ptope.nodes[1] = 1;
        ptope.nodes[2] = 2;
        ptope.nodes[3] = 3;
        ptope.nodes[4] = 4;

        std::vector<Vector4r> points(
            {{0., 0., 0., 0.}, {2., 0., 0., 0.}, {0., 2., 0., 0.}, {0., 0., 2., 0.}, {0., 0., 0., 2.}});

        SimplexInterpolator<4> interp;
        std::array<Pentatope4, 16> sub_simplices;
        std::vector<Vector4r> sub_points;

        red_refinement<4, 4, 16>(ptope, points, sub_simplices, sub_points, interp);
    }

    template <mars::Integer Dim, mars::Integer ManifoldDim>
    void test_mesh(mars::Mesh<Dim, ManifoldDim> &mesh) {
        using namespace mars;
        using MeshD = mars::Mesh<Dim, ManifoldDim>;

        RedGreenRefinement<MeshD> rgr(mesh);

        static const Integer n_refinements = 1;
        mesh.build_dual_graph();
        Integer nbs = mesh.n_boundary_sides();
        std::cout << "n_boundary_sides: " << nbs << std::endl;
        mesh.check_side_ordering();
        // mesh.describe_dual_graph(std::cout);

        std::cout << "-------------------------" << std::endl;

        // mesh.red_refine_element(0);
        rgr.uniformly_refine(n_refinements);
        std::cout << "n_elements: " << mesh.n_active_elements() << std::endl;
        std::cout << "n_nodes:    " << mesh.n_nodes() << std::endl;

        mesh.build_dual_graph();
        mesh.check_side_ordering();
        std::cout << "n_boundary_sides: " << mesh.n_boundary_sides()
                  << " == " << Power<2, (ManifoldDim - 1) * n_refinements>::value * nbs << std::endl;

        // mesh.describe(std::cout, true);
        // mesh.describe(std::cout, false);

        VTKMeshWriter<MeshD> w;
        w.write("mesh_red_refined_benchmark.vtu", mesh);
    }

    template <mars::Integer Dim, mars::Integer ManifoldDim>
    void generate_mesh(mars::Mesh<Dim, ManifoldDim> &mesh, int level) {
        using namespace mars;
        using MeshD = mars::Mesh<Dim, ManifoldDim>;

        RedGreenRefinement<MeshD> rgr(mesh);

        static const Integer n_refinements = 1;
        mesh.build_dual_graph();
        Integer nbs = mesh.n_boundary_sides();
        std::cout << "n_boundary_sides: " << nbs << std::endl;
        mesh.check_side_ordering();
        // mesh.describe_dual_graph(std::cout);

        std::cout << "-------------------------" << std::endl;

        // mesh.red_refine_element(0);
        for (int i = 0; i < level; i++) rgr.uniformly_refine(n_refinements);
        std::cout << "n_elements: " << mesh.n_active_elements() << std::endl;
        std::cout << "n_nodes:    " << mesh.n_nodes() << std::endl;

        mesh.build_dual_graph();
        mesh.check_side_ordering();
        std::cout << "n_boundary_sides: " << mesh.n_boundary_sides()
                  << " == " << Power<2, (ManifoldDim - 1) * n_refinements>::value * nbs << std::endl;

        // mesh.describe(std::cout, true);
        // mesh.describe(std::cout, false);

        VTKMeshWriter<MeshD> w;
        w.write("mesh_red_refined_uniformtest" + std::to_string(mesh.Dim) + ".vtu", mesh);
    }

    template <mars::Integer Dim, mars::Integer ManifoldDim>
    static void test_generate_3D(mars::Mesh<Dim, ManifoldDim> &mesh) {
        using namespace mars;

        std::cout << "======================================\n";
        using MeshD = mars::Mesh<Dim, ManifoldDim>;

        RedGreenRefinement<MeshD> rgr(mesh);
        rgr.red_refine({0});
        // mesh.describe(std::cout);
        write_element("elem_3.eps", rgr, 0, 10, INVALID_INDEX);
        write_element_with_sides("elem_sides_3.eps", mesh, 0, 10, INVALID_INDEX);

        write_element_with_subsurfaces("elem_ss_3.eps", mesh, 0, 10);
        std::cout << "======================================\n";

        VTKMeshWriter<MeshD> w;
        w.write("mesh_red_refined_test" + std::to_string(mesh.Dim) + ".vtu", mesh);
    }

    static void test_mfem_mesh_2D() {
        using namespace mars;
        std::cout << "======================================\n";
        Mesh2 mesh;
        read_mesh("../data/square_2.MFEM", mesh);

        RedGreenRefinement<Mesh2> rgr(mesh);

        test_mesh(mesh);
        write_mesh("mesh_2.eps", mesh, 10., PLOT_ID);

        std::cout << "mesh" << std::endl;
        mesh.describe_boundary_elements(std::cout);

        rgr.red_refine({2, 4});

        std::cout << "red 1" << std::endl;
        mesh.describe_boundary_elements(std::cout);
        write_mesh("mesh_2_r1.eps", mesh, 10., PLOT_PARENT_TAG);

        rgr.green_refine();

        std::cout << "green 1" << std::endl;
        mesh.describe_boundary_elements(std::cout);
        write_mesh("mesh_2_rg1.eps", mesh, 10., PLOT_PARENT_TAG);

        write_element("elem_2.eps", rgr, 2, 10, INVALID_INDEX);
        write_element_with_sides("elem_sides_2.eps", mesh, 0, 10, INVALID_INDEX);

        write_element_with_subsurfaces("elem_ss_2.eps", mesh, 2, 10);
    }

    static void test_mfem_mesh_3D() {
        using namespace mars;

        std::cout << "======================================\n";
        Mesh3 mesh;
        read_mesh("../data/cube_6.MFEM", mesh, true);
        // test_mesh(mesh);
        RedGreenRefinement<Mesh3> rgr(mesh);
        rgr.red_refine({0, 2});
        // mesh.describe(std::cout);
        write_element("elem_3.eps", rgr, 0, 10, INVALID_INDEX);
        write_element_with_sides("elem_sides_3.eps", mesh, 0, 10, INVALID_INDEX);

        write_element_with_subsurfaces("elem_ss_3.eps", mesh, 0, 10);
        std::cout << "======================================\n";

        VTKMeshWriter<Mesh3> w;
        w.write("mesh_red_refined_benchmark.vtu", mesh);
    }

    static void test_mfem_mesh_4D() {
        using namespace mars;
        std::cout << "======================================\n";
        Mesh4 mesh;
        read_mesh("../data/pentatope_1.MFEM", mesh);
        // read_mesh("../data/cube4d_24.MFEM", mesh);
        // test_mesh(mesh);

        RedGreenRefinement<Mesh4> rgr(mesh);
        rgr.red_refine({0});
        // mesh.describe(std::cout);
        // mesh.describe_dual_graph(std::cout);
        write_element("elem_4.eps", rgr, 0, 10, INVALID_INDEX);
        write_element_with_sides("elem_sides_4.eps", mesh, 6, 10, INVALID_INDEX);
        std::cout << "======================================\n";

        MultilevelElementMap<3, 2> mlem;
        mlem.update(mesh.elem(0));
        mlem.update(mesh.elem(1));
        // mlem.describe(std::cout);

        write_element_with_subsurfaces("elem_ss_4.eps", mesh, 0, 10);
    }

    static void test_primitives() {
        Tetrahedron4 tet;
        tet.nodes[0] = 0;
        tet.nodes[1] = 1;
        tet.nodes[2] = 2;
        tet.nodes[3] = 3;

        Triangle4 tri4;
        tri4.nodes[0] = 1;
        tri4.nodes[1] = 2;
        tri4.nodes[2] = 3;

        Line2 line;
        line.nodes[0] = 0;
        line.nodes[1] = 1;

        Triangle2 tri2;
        tri2.nodes[0] = 0;
        tri2.nodes[1] = 1;
        tri2.nodes[2] = 2;

        Pentatope4 ptope;
        ptope.nodes[0] = 0;
        ptope.nodes[1] = 1;
        ptope.nodes[2] = 2;
        ptope.nodes[3] = 3;
        ptope.nodes[4] = 4;

        std::vector<Vector4r> points4(
            {{0., 0., 0., 0.}, {2., 0., 0., 0.}, {0., 2., 0., 0.}, {0., 0., 2., 0.}, {0., 0., 0., 2.}});

        std::vector<Vector2r> points2({{0., 0.}, {1., 0.}, {0., 1.}});

        std::cout << "===================" << std::endl;
        std::cout << "vol(line2): " << volume(line, points2) << std::endl;
        std::cout << "vol(tri2):  " << volume(tri2, points2) << std::endl;
        std::cout << "vol(tri4):  " << volume(tri4, points4) << std::endl;
        std::cout << "vol(tet4):  " << volume(tet, points4) << std::endl;
        std::cout << "vol(ptope): " << volume(ptope, points4) << std::endl;
        std::cout << "===================" << std::endl;

        Matrix<Real, 4, 4> J44;
        jacobian(ptope, points4, J44);
        std::cout << "J(ptope):\n" << J44 << std::endl;

        Matrix<Real, 4, 2> Jtri4;
        jacobian(tri4, points4, Jtri4);
        check_and_fix_jac(Jtri4);
        std::cout << "J(tri4):\n" << Jtri4 << std::endl;

        auto n1 = normal(line, points2);
        auto n3 = normal(tri4, points4);
        auto n4 = normal(tet, points4);

        std::cout << n1 << std::endl;
        std::cout << n3 << std::endl;
        std::cout << n4 << std::endl;
    }

    void test_det() {
        static const Integer N = 6;
        Matrix<Real, N, N> mat;

        for (Integer i = 0; i < N; ++i) {
            mat(i, i) = 1.;
        }

        Real d = det(mat);
        assert(abs(d - 1.) < 1e-16);

        //////////////////////////////////////////////////////

        Matrix<Real, 5, 5> magic5(
            {17, 24, 1, 8, 15, 23, 5, 7, 14, 16, 4, 6, 13, 20, 22, 10, 12, 19, 21, 3, 11, 18, 25, 2, 9});

        Real det_magic5 = det(magic5);
        assert(abs(det_magic5 - 5069999.999999999) < 1e-8);

        //////////////////////////////////////////////////////

        Matrix<Real, 6, 6> magic6({35 + 1, 1, 6,     26, 19,     24, 3, 32 + 1, 7,  21,     23, 25,
                                   31,     9, 2 + 1, 22, 27,     20, 8, 28,     33, 17 + 1, 10, 15,
                                   30,     5, 34,    12, 14 + 1, 16, 4, 36,     29, 13,     18, 11 + 1});

        Real det_magic6 = det(magic6);
        assert(abs(det_magic6 - 7745920.00000) < 1e-6);

        //////////////////////////////////////////////////////
    }

    void stream_write_read() {
        std::ostringstream oss;
        std::istringstream iss;

        Integer num = 10;
        Integer read_num = 0;

        write(num, oss);
        iss.str(oss.str());
        read(read_num, iss);

        assert(num == read_num);
    }

    void run_tests(const int level) {
        // test_det();
        //  test_midpoint_index();
        //  test_red_refinement_interpolator();
        //  test_red_refinement();
        // test_mfem_mesh_3D();
        //  test_mfem_mesh_4D();
        //  test_mfem_mesh_2D();
        /*Mesh3 mesh;
        read_mesh("../data/write/tetrakis.MFEM", mesh, true);
        generate_mesh<3,3>(mesh,1);
        VTKMeshWriter<Mesh3> w;
                                w.write("cube_red"+ std::to_string(mesh.Dim) +".vtu", mesh);*/

        Mesh<3, 2> tri;
        read_mesh("../data/square_3.MFEM", tri, true);
        generate_mesh(tri, level);

        // test_generate_3D(tri);

        // write_mesh_MFEM("../data/write/cube_6.MFEM",mesh);

        // test_mfem_mesh_2D();
    }
}  // namespace mars
