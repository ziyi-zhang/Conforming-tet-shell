#include <shell/common.hpp>
#include <shell/Label.h>

#include <igl/copyleft/cgal/orient3D.h>
#include <Eigen/Dense>

#include <array>


namespace tetshell {

namespace {

bool point_in_tetrahedron(const Vec3d& point, const Vec3d& T0, const Vec3d& T1, const Vec3d& T2, const Vec3d& T3) {

    using igl::copyleft::cgal::orient3D;
    return  orient3D(T0.data(), T3.data(), T1.data(), point.data()) >= 0 &&
            orient3D(T1.data(), T3.data(), T2.data(), point.data()) >= 0 &&
            orient3D(T0.data(), T1.data(), T2.data(), point.data()) >= 0 &&
            orient3D(T0.data(), T2.data(), T3.data(), point.data()) >= 0;
}


bool point_in_prism(const Vec3d& point, bool tetra_split_AB, const std::array<Vec3d, 6>& verts) {

    auto tets = tetra_split_AB ? TETRA_SPLIT_A : TETRA_SPLIT_B;

    for (int i = 0; i < 3; i++)
        if (point_in_tetrahedron(point, verts[tets[i][0]], verts[tets[i][1]],
                                verts[tets[i][2]], verts[tets[i][3]]))
        return true;
    return false;
}


Vec3d GetBarycenter(const Eigen::MatrixXd &VO, const Eigen::MatrixXi &TO, int i) {

    Vec3d center;

    center(0) = (VO(TO(i, 0), 0) + VO(TO(i, 1), 0) + VO(TO(i, 2), 0) + VO(TO(i, 3), 0)) / 4.0;
    center(1) = (VO(TO(i, 0), 1) + VO(TO(i, 1), 1) + VO(TO(i, 2), 1) + VO(TO(i, 3), 1)) / 4.0;
    center(2) = (VO(TO(i, 0), 2) + VO(TO(i, 1), 2) + VO(TO(i, 2), 2) + VO(TO(i, 3), 2)) / 4.0;
    return center;
}


bool PointInShell(const Vec3d &center, const shell_t &shell) {

    int N = shell.size();

    for (int i=0; i<N; i++) {
        if (point_in_prism(center, true, shell[i]))
            return true;
    }
    return false;
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LabelTet(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI, const Eigen::MatrixXd &VO, const Eigen::MatrixXi &TO, Eigen::VectorXi &labels) {

    const int M = VI.rows() / 4;
    const int Ntri = FI.rows() / 4;
    const int Ntet = TO.rows();
    labels.setZero(Ntet, 1);

    // Get four surfaces from VI (temp)
    Eigen::MatrixXd V_inner, V_bottom, V_top, V_outer;
    Eigen::MatrixXi F;
    V_inner = VI.block(0, 0, M, 3);
    V_bottom = VI.block(M, 0, M, 3);
    V_top = VI.block(M*2, 0, M, 3);
    V_outer = VI.block(M*3, 0, M, 3);
    F = FI.block(0, 0, Ntri, 3);
    // Retrieve prisms from dual-shell
    shell_t shell_inner_bottom, shell_bottom_top, shell_top_outer;
    ConstructShellFromTri(V_inner, V_bottom, F, shell_inner_bottom);
    ConstructShellFromTri(V_bottom, V_top, F, shell_bottom_top);
    ConstructShellFromTri(V_top, V_outer, F, shell_top_outer);

    // now loop over all tets and label them
    // dumb way, not optimized
    for (int i=0; i<Ntet; i++) {
        
        Vec3d center = GetBarycenter(VO, TO, i);
        if (PointInShell(center, shell_inner_bottom))
            labels(i) = 1;
        else if (PointInShell(center, shell_bottom_top))
            labels(i) = 2;
        else if (PointInShell(center, shell_top_outer))
            labels(i) = 3;
    }
}


void ConstructShellFromTri(const Eigen::MatrixXd &V_bottom, const Eigen::MatrixXd &V_top, const Eigen::MatrixXi &F, shell_t &shell) {

    const int N = F.rows();

    shell.clear();
    for (int i=0; i<N; i++) {

        prism_t prism;
        prism[0] = Vec3d(V_bottom(F(i, 0), 0), V_bottom(F(i, 0), 1), V_bottom(F(i, 0), 2));
        prism[1] = Vec3d(V_bottom(F(i, 1), 0), V_bottom(F(i, 1), 1), V_bottom(F(i, 1), 2));
        prism[2] = Vec3d(V_bottom(F(i, 2), 0), V_bottom(F(i, 2), 1), V_bottom(F(i, 2), 2));
        prism[3] = Vec3d(V_top(F(i, 0), 0), V_top(F(i, 0), 1), V_top(F(i, 0), 2));
        prism[4] = Vec3d(V_top(F(i, 1), 0), V_top(F(i, 1), 1), V_top(F(i, 1), 2));
        prism[5] = Vec3d(V_top(F(i, 2), 0), V_top(F(i, 2), 1), V_top(F(i, 2), 2));
        shell.push_back(prism);
    }
}


}  // namespace tetshell
