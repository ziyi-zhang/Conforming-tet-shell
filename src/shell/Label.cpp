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
#include <algorithm>
#include <string>


namespace tetshell {

using namespace std;
using tetwild::Point_3;
using tetwild::Segment_3;
using tetwild::logger;

namespace {

Point_3 GetBarycenter(const std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO, int i) {

    Point_3 center;

    center = CGAL::centroid(VO[TO[i][0]].pos, VO[TO[i][1]].pos, VO[TO[i][2]].pos, VO[TO[i][3]].pos);
    /*
    // DEBUG PURPOSE
        double pos_0_x = CGAL::to_double(VO[TO[i][0]].pos.x());
        double pos_0_y = CGAL::to_double(VO[TO[i][0]].pos.y());
        double pos_0_z = CGAL::to_double(VO[TO[i][0]].pos.z());

        double pos_1_x = CGAL::to_double(VO[TO[i][1]].pos.x());
        double pos_1_y = CGAL::to_double(VO[TO[i][1]].pos.y());
        double pos_1_z = CGAL::to_double(VO[TO[i][1]].pos.z());

        double pos_2_x = CGAL::to_double(VO[TO[i][2]].pos.x());
        double pos_2_y = CGAL::to_double(VO[TO[i][2]].pos.y());
        double pos_2_z = CGAL::to_double(VO[TO[i][2]].pos.z());

        double pos_3_x = CGAL::to_double(VO[TO[i][3]].pos.x());
        double pos_3_y = CGAL::to_double(VO[TO[i][3]].pos.y());
        double pos_3_z = CGAL::to_double(VO[TO[i][3]].pos.z());

        double center_x = CGAL::to_double(center.x());
        double center_y = CGAL::to_double(center.y());
        double center_z = CGAL::to_double(center.z());
    */
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
    /*
    // DEBUG PURPOSE
    double T0_x = CGAL::to_double(T0.x());
    double T0_y = CGAL::to_double(T0.y());
    double T0_z = CGAL::to_double(T0.z());
    double T1_x = CGAL::to_double(T1.x());
    double T1_y = CGAL::to_double(T1.y());
    double T1_z = CGAL::to_double(T1.z());
    double T2_x = CGAL::to_double(T2.x());
    double T2_y = CGAL::to_double(T2.y());
    double T2_z = CGAL::to_double(T2.z());
    double T3_x = CGAL::to_double(T3.x());
    double T3_y = CGAL::to_double(T3.y());
    double T3_z = CGAL::to_double(T3.z());
    */

    CGAL::Tetrahedron_3<tetwild::K> tet(T0, T1, T2, T3);
    if (tet.is_degenerate()) return false;  // in case of singularity
    CGAL::Oriented_side side = tet.oriented_side(point);
    CGAL::Orientation ori = CGAL::orientation(T0, T1, T2, T3);
    if (ori == CGAL::POSITIVE)
        return (side == CGAL::Oriented_side::ON_POSITIVE_SIDE) || (side == CGAL::ON_ORIENTED_BOUNDARY);
    else {
        // tetwild::log_and_throw("point_in_tetrahedron: flipped tet detected.");
        return (side == CGAL::Oriented_side::ON_NEGATIVE_SIDE) || (side == CGAL::ON_ORIENTED_BOUNDARY);
    }
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
        if (point_in_prism(center, shell[i][1] > shell[i][2], prism_vertices)) {
            // if (debugFlag) std::cerr << "debugFlag " << i << std::endl;
            return true;
        }
    }
    return false;
}


int TetRegion(
    const std::vector<tetwild::TetVertex> &VO, 
    const std::vector<std::array<int, 4>> &TO,
    const std::vector<Point_3> &VI, 
    const DualShell_t &dualShell,
    const int tetIdx) {

    Point_3 center = GetBarycenter(VO, TO, tetIdx);
    if (PointInShell(center, dualShell.shell_inner_bottom, VI))
        return SHELL_INNER_BOTTOM;
    else if (PointInShell(center, dualShell.shell_bottom_top, VI))
        return SHELL_BOTTOM_TOP;
    else if (PointInShell(center, dualShell.shell_top_outer, VI))
        return SHELL_TOP_OUTER;
    
    return 0;  // default zero, not in any layer of the shell
}


int GetRegionType(int oldRegionType, int surfaceType) {

    // do not change region type unless met a shell surface
    // The shell must be closed
    if (surfaceType == NOT_SUR) return oldRegionType;

    const int NON_SHELL_REGION = 0;
    if (oldRegionType == NON_SHELL_REGION) {
        switch (surfaceType) {
            case SURFACE_INNER:
                return SHELL_INNER_BOTTOM;
            case SURFACE_OUTER:
                return SHELL_TOP_OUTER;
            default:
                std::cerr << surfaceType << std::endl;
                tetwild::log_and_throw("GetRegionType: exception in NON_SHELL_REGION.");
                break;
        }
    } else if (oldRegionType == SHELL_INNER_BOTTOM) {
        switch (surfaceType) {
            case SURFACE_INNER:
                return NON_SHELL_REGION;
            case SURFACE_BOTTOM:
                return SHELL_BOTTOM_TOP;
            default:
                std::cerr << surfaceType << std::endl;
                tetwild::log_and_throw("GetRegionType: exception in SHELL_INNER_BOTTOM.");
                break;
        }
    } else if (oldRegionType == SHELL_BOTTOM_TOP) {
        switch (surfaceType) {
            case SURFACE_BOTTOM:
                return SHELL_INNER_BOTTOM;
            case SURFACE_TOP:
                return SHELL_TOP_OUTER;
            default:
                std::cerr << surfaceType << std::endl;
                tetwild::log_and_throw("GetRegionType: exception in SHELL_BOTTOM_TOP.");
                break;
        }
    } else if (oldRegionType == SHELL_TOP_OUTER) {
        switch (surfaceType) {
            case SURFACE_TOP:
                return SHELL_BOTTOM_TOP;
            case SURFACE_OUTER:
                return NON_SHELL_REGION;
            default:
                std::cerr << surfaceType << std::endl;
                tetwild::log_and_throw("GetRegionType: exception in SHELL_TOP_OUTER.");
                break;
        }
    }

    tetwild::log_and_throw("GetRegionType: invalid input.");
    return 0;  // never reach here
}


void LabelTet(
    const tetwild::Args &args,
    const std::vector<Point_3> &VI, 
    const Eigen::MatrixXi &FI, 
    const std::vector<tetwild::TetVertex> &VO, 
    const std::vector<std::array<int, 4>> &TO,
    const std::vector<std::array<int, 4>> &face_on_shell, 
    DualShell_t &dualShell,
    Eigen::VectorXi &labels) {

    const int Ntet = TO.size();
    labels.setOnes(Ntet, 1);  
    labels = labels * (-1);  // default -1 (unlabeled). 0 means not in any layer of the shell.

    // Get four surfaces from VI
    GenDualShell(VI, FI, dualShell);

    if (args.brute_label) {
    // brute label

        // loop over all tets and check which prism contains them
        for (int i=0; i<Ntet; i++) {
            labels[i] = TetRegion(VO, TO, VI, dualShell, i);
        }
    } else {
    // BFS label

        // BFS to label tets
        std::queue<std::pair<int, int> > Q;  // <tetIdx, regionType>
                                            // Q stores tets that have been visited & labeled but can be used to expand search
        std::unordered_set<int> uset, set_tmp;
        std::vector<bool> visited(Ntet, false);  // whether this tet has been visited and labeled

        // first element
        int regionTet0 = TetRegion(VO, TO, VI, dualShell, 0);
        Q.push(std::make_pair(0, regionTet0));
        labels[0] = regionTet0;
        visited[0] = true;

        // Start BFS
        int cnt = 0;  // brute_label_validation only
        while (!Q.empty()) {

            // tetIdx-th tetrahedron
            int oldTetIdx = Q.front().first;
            int oldRegionType = Q.front().second;
            Q.pop();

            // Consider four faces seperately
            for (int i=0; i<4; i++) {

                int vert1 = TO[oldTetIdx][(0+i) % 4];
                int vert2 = TO[oldTetIdx][(1+i) % 4];
                int vert3 = TO[oldTetIdx][(2+i) % 4];
                int oppositeVert = (3+i) % 4;  // not the index, only need to know which one among the four
                // Find the intersection of {vert1.conn_tets, vert2.conn_tets, vert3.conn_tets}
                UnorderedsetIntersection(VO[vert1].conn_tets, VO[vert2].conn_tets, set_tmp);
                UnorderedsetIntersection(set_tmp, VO[vert3].conn_tets, uset);

                for (int newTetIdx : uset) {

                    if (newTetIdx == oldTetIdx) continue;
                    if (visited[newTetIdx]) continue;
                    // decide the label of newTetIdx
                    int surfaceType = face_on_shell[oldTetIdx][oppositeVert];
                    int newRegionType = GetRegionType(oldRegionType, surfaceType);
                    labels[newTetIdx] = newRegionType;
                    visited[newTetIdx] = true;

                    // DEBUG ONLY
                    if (args.brute_label_validation) {
                        if (newRegionType != TetRegion(VO, TO, VI, dualShell, newTetIdx)) {
                            labels[newTetIdx] = 5;  // err code
                            cnt++;
                            logger().error("Wrong tet label: newRegionType = {}, tetRegion = {}", newRegionType, TetRegion(VO, TO, VI, dualShell, newTetIdx));
                            std::string errCode = "newTetIdx = " + std::to_string(newTetIdx) + "  oldRegionType = " + std::to_string(oldRegionType) + "  surfaceType = " + std::to_string(surfaceType);
                            std::cout << errCode << std::endl;
                            // return;
                        }
                    }

                    // push to Q for further search
                    Q.push(std::make_pair(newTetIdx, newRegionType));
                }
            }
        }

        // assert
        /// Note: even the model is not connected, the use of a bounding box will make them connected
        if (std::count(visited.begin(), visited.end(), false))
            tetwild::log_and_throw("LabelTet: not all tets are visited.");
        if (args.brute_label_validation)  // DEBUG PURPOSE
            std::cerr << "brute_label_validation mismatch cnt = " << cnt << std::endl;
    }

/*
    // this code can be used to visualize "face_on_shell"
    for (int i=0; i<Ntet; i++) {

        const std::array<int, 4> &faceInfo = face_on_shell[i];
        int countNonZero = int(faceInfo[0] > 0) + int(faceInfo[1] > 0) + int(faceInfo[2] > 0) + int(faceInfo[3] > 0);
        if (countNonZero == 0)
            labels[i] = -1;
        else if (countNonZero > 1)
            labels[i] = 5;
        else if (faceInfo[0] > 0)
            labels[i] = faceInfo[0];
        else if (faceInfo[1] > 0)
            labels[i] = faceInfo[1];
        else if (faceInfo[2] > 0)
            labels[i] = faceInfo[2];
        else if (faceInfo[3] > 0)
            labels[i] = faceInfo[3];
    }
*/
    logger().info("Tetrahedra label done.");
}

}  // namespace tetshell
