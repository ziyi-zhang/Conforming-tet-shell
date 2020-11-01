#include <shell/common.hpp>
#include <shell/Label.h>
#include <shell/Shell.h>
#include <shell/Utils.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>
#include <tetwild/TetmeshElements.h>

#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <igl/copyleft/cgal/orient3D.h>
#include <Eigen/Dense>
#include <map>
#include <set>
#include <unordered_set>
#include <array>
#include <vector>
#include <queue>


namespace tetshell {

using namespace std;
using tetwild::Point_3;
using tetwild::Segment_3;
using tetwild::logger;

namespace {

Point_3 GetBarycenter(const std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO, int i) {

    Point_3 center;

    center = CGAL::centroid(VO[TO[i][0]].pos, VO[TO[i][1]].pos, VO[TO[i][2]].pos, VO[TO[i][3]].pos);
    return center;
}


/// Concatenate TetMesh 2 under TetMesh 1
/// This may create duplicate vertices
/// [ DEPRECATED ]
void UnionTetMesh(
    std::vector<tetwild::TetVertex> &V1, 
    std::vector<std::array<int, 4>> &T1, 
    const std::vector<tetwild::TetVertex> &V2, 
    const std::vector<std::array<int, 4>> &T2) {

    const int N = V1.size();

    // concatenate V1 V2 to V1
    V1.insert( V1.end(), V2.begin(), V2.end() );

    // T = T2 + N
    std::vector<std::array<int, 4>> T = T2;
    for (auto it=T.begin(); it!=T.end(); it++) {
        it->at(0) += N;
        it->at(1) += N;
        it->at(2) += N;
        it->at(3) += N;
    }

    // concatenate T1 T to T1
    T1.insert( T1.end(), T.begin(), T.end() );
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////

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


bool PointInShell(const Point_3 &center, const shell_t &shell, const std::vector<Point_3> &VI) {

    int N = shell.size();

    for (int i=0; i<N; i++) {

        std::array<Point_3, 6> prism_vertices;
        prism_vertices[0] = VI[shell[i][0]];
        prism_vertices[1] = VI[shell[i][1]];
        prism_vertices[2] = VI[shell[i][2]];
        prism_vertices[3] = VI[shell[i][3]];
        prism_vertices[4] = VI[shell[i][4]];
        prism_vertices[5] = VI[shell[i][5]];
        if (point_in_prism(center, true, prism_vertices))  // fine to always use TYPE-A
            return true;
    }
    return false;
}


void LabelTet(
    const std::vector<Point_3> &VI, 
    const Eigen::MatrixXi &FI, 
    const std::vector<tetwild::TetVertex> &VO, 
    const std::vector<std::array<int, 4>> &TO,
    const std::vector<std::array<int, 4>> &face_on_shell, 
    DualShell_t &dualShell,
    Eigen::VectorXi &labels) {

    const int Ntet = TO.size();
    labels.setZero(Ntet, 1);  // default zero, not in any layer of the shell

    // Get four surfaces from VI
    GenDualShell(VI, FI, dualShell);

    // loop over all tets and check which prism contains them
    /*
    for (int i=0; i<Ntet; i++) {

        Point_3 center = GetBarycenter(VO, TO, i);
        if (PointInShell(center, dualShell.shell_inner_bottom, VI))
            labels(i) = 1;
        else if (PointInShell(center, dualShell.shell_bottom_top, VI))
            labels(i) = 2;
        else if (PointInShell(center, dualShell.shell_top_outer, VI))
            labels(i) = 3;
    }
    */

    // BFS to label tets
    std::queue<int> Q;
    

    logger().info("Tetrahedra label done.");
}

}  // namespace tetshell
