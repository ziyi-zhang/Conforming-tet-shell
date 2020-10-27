#include <shell/common.hpp>
#include <shell/Shell.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>

#include <Eigen/Dense>


namespace tetshell {

using tetwild::Point_3;

namespace {

void ConstructShellFromTri(const std::vector<Point_3> &V_bottom, const std::vector<Point_3> &V_top, const RowMatX3i &F, shell_t &shell) {

    const int N = F.rows();

    shell.clear();
    for (int i=0; i<N; i++) {

        prism_t prism;
        prism[0] = V_bottom[F(i, 0)];
        prism[1] = V_bottom[F(i, 1)];
        prism[2] = V_bottom[F(i, 2)];
        prism[3] = V_top   [F(i, 0)];
        prism[4] = V_top   [F(i, 1)];
        prism[5] = V_top   [F(i, 2)];
        shell.push_back(prism);
    }
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////

void GenDualShell(const std::vector<Point_3> &VI, const Eigen::MatrixXi &FI, DualShell_t &dualShell) {

    const int M = VI.size() / 4;
    const int Ntri = FI.rows() / 4;

    // assert size - check this even in release
    if ((VI.size() % 4 != 0) ||
        (FI.rows() % 4 != 0)) {
        assert(false);
        std::cerr << "Wrong input size" << std::endl;
        tetwild::log_and_throw("GenDualShell: Wrong input size");
    }
    // assert F - check this even in release
    Eigen::MatrixXi F1, F2, F3, F4;
    F1 = FI.block(0, 0, Ntri, 3);
    F2 = FI.block(Ntri, 0, Ntri, 3);
    F3 = FI.block(2*Ntri, 0, Ntri, 3);
    F4 = FI.block(3*Ntri, 0, Ntri, 3);
    if ((! F1.isApprox(F2.array() - M)) ||
        (! F1.isApprox(F3.array() - M*2)) ||
        (! F1.isApprox(F4.array() - M*3))) {
        assert(false);
        std::cerr << "Wrong input index matrix" << std::endl;
        tetwild::log_and_throw("GenDualShell: Wrong input index matrix");
    }

    // split input
    dualShell.F = F1;
    dualShell.V_inner = {VI.begin(), VI.begin()+M};
    dualShell.V_bottom = {VI.begin()+M, VI.begin()+2*M};
    dualShell.V_top = {VI.begin()+2*M, VI.begin()+3*M};
    dualShell.V_outer = {VI.begin()+3*M, VI.begin()+4*M};

    // Retrieve prisms from dual-shell
    ConstructShellFromTri(dualShell.V_inner, dualShell.V_bottom, F1, dualShell.shell_inner_bottom);
    ConstructShellFromTri(dualShell.V_bottom, dualShell.V_top, F1, dualShell.shell_bottom_top);
    ConstructShellFromTri(dualShell.V_top, dualShell.V_outer, F1, dualShell.shell_top_outer);
}


bool point_in_tetrahedron(const Point_3& point, const Point_3& T0, const Point_3& T1, const Point_3& T2, const Point_3& T3) {

    /*
    // using libigl::orient3d
    using igl::copyleft::cgal::orient3D;
    return  orient3D(T0.data(), T3.data(), T1.data(), point.data()) >= 0 &&
            orient3D(T1.data(), T3.data(), T2.data(), point.data()) >= 0 &&
            orient3D(T0.data(), T1.data(), T2.data(), point.data()) >= 0 &&
            orient3D(T0.data(), T2.data(), T3.data(), point.data()) >= 0;
    */

    CGAL::Tetrahedron_3<tetwild::K> tet(T0, T1, T2, T3);
    CGAL::Oriented_side side = tet.oriented_side(point);
    CGAL::Orientation ori = CGAL::orientation(T0, T1, T2, T3);
    if (ori == CGAL::POSITIVE)
        return side == CGAL::Oriented_side::ON_POSITIVE_SIDE;
    else
        return side == CGAL::Oriented_side::ON_NEGATIVE_SIDE;
}


bool point_in_prism(const Point_3& point, bool tetra_split_AB, const std::array<Point_3, 6>& verts) {

    auto tets = tetra_split_AB ? TETRA_SPLIT_A : TETRA_SPLIT_B;

    for (int i = 0; i < 3; i++)
        if (point_in_tetrahedron(point, verts[tets[i][0]], verts[tets[i][1]],
                                verts[tets[i][2]], verts[tets[i][3]]))
        return true;
    return false;
} 


bool PointInShell(const Point_3 &center, const shell_t &shell) {

    int N = shell.size();

    for (int i=0; i<N; i++) {
        if (point_in_prism(center, true, shell[i]))  // fine to always use TYPE-A
            return true;
    }
    return false;
}

}  // namespace tetshell
