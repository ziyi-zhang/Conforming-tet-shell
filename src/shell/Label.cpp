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
#include <array>
#include <vector>


namespace tetshell {

using namespace std;
using tetwild::Point_3;
using tetwild::logger;

namespace {

Point_3 GetBarycenter(const std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO, int i) {

    Point_3 center;

    center = CGAL::centroid(VO[TO[i][0]].pos, VO[TO[i][1]].pos, VO[TO[i][2]].pos, VO[TO[i][3]].pos);
    return center;
}


/// Concatenate TetMesh 2 under TetMesh 1
/// This may create duplicate vertices
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


/// Map index of VI to index of VO for all indices in "shell"
void MapIndex(
    const std::vector<tetwild::TetVertex> &VO,  // VO
    const DualShell_t &dualShell,  // VI
    const shell_t &shell, 
    std::map<int, int> &idxMap) {

    const int N = shell.size();
    const std::vector<tetwild::Point_3> &VI = dualShell.V;
    idxMap.clear();

    for (int i=0; i<N; i++) {

        const prism_t &prism = shell[i];
        for (int j=0; j<6; j++) {

            int VI_index = prism[j];
            if (idxMap.find(VI_index) == idxMap.end()) {
                // no mapping exist yet
                bool found = false;
                for (auto it=VO.begin(); it!=VO.end(); it++) {
                    if (it->pos == VI[VI_index]) {
                        idxMap.insert(std::make_pair(VI_index, it-VO.begin()));
                        found = true;
                        break;
                    }
                }
                if (found) continue;
                // Should not reach here
                tetwild::log_and_throw("MapIndex: A vertex in VI not found in VO");
            }
        }
    }
}


bool IsPointInsideTri(
    const Point_3 &pt1, 
    const Point_3 &base_pt1, 
    const Point_3 &base_pt2, 
    const Point_3 &base_pt3) {

    tetwild::Triangle_3 tri(base_pt1, base_pt2, base_pt3);
    return tri.has_on(pt1);
}


/// Is triangle (pt1, pt2, pt3) inside triangle (base_pt1, base_pt2, base_pt3)
bool IsTriInsideTri(
    const Point_3 &pt1, 
    const Point_3 &pt2, 
    const Point_3 &pt3, 
    const Point_3 &base_pt1, 
    const Point_3 &base_pt2, 
    const Point_3 &base_pt3) {

    bool pt1_good = IsPointInsideTri(pt1, base_pt1, base_pt2, base_pt3);
    bool pt2_good = IsPointInsideTri(pt2, base_pt1, base_pt2, base_pt3);
    bool pt3_good = IsPointInsideTri(pt3, base_pt1, base_pt2, base_pt3);

    return pt1_good && pt2_good && pt3_good;
}


void CleanTetMesh(
    const std::vector<bool> &t_is_removed,
    std::vector<tetwild::TetVertex> &V, 
    std::vector<std::array<int, 4>> &T, 
    Eigen::VectorXi &labels, 
    std::vector<std::array<int, 4>> &is_surface_facet, 
    std::vector<std::array<int, 4>> &face_on_shell) {

    std::vector<tetwild::TetVertex> V_out;
    std::vector<std::array<int, 4>> T_out;
    Eigen::VectorXi labels_out;
    std::vector<std::array<int, 4>> is_surface_facet_out;
    std::vector<std::array<int, 4>> face_on_shell_out;
    const int tetNum = std::count(t_is_removed.begin(), t_is_removed.end(), false);

    // v_ids is the vector of index of vertices
    std::vector<int> v_ids;
    for (int i=0; i<T.size(); i++) {
        if (t_is_removed[i])
            continue;
        for (int j=0; j<4; j++)
            v_ids.push_back(T[i][j]);
    }
    // unique "v_ids"
    std::sort(v_ids.begin(), v_ids.end());
    v_ids.erase(std::unique(v_ids.begin(), v_ids.end()), v_ids.end());
    // "map_ids" is defined as this for all used vertices:
    // Key: index in V
    // Value: new index i
    std::unordered_map<int, int> map_ids;
    for (int i = 0; i < v_ids.size(); i++)
        map_ids[v_ids[i]] = i;

    // Prepare 
    V_out.clear();
    T_out.clear();
    labels_out.resize(tetNum, 1);
    is_surface_facet_out.clear();
    face_on_shell_out.clear();
    // Fill V
    for (int i = 0; i < v_ids.size(); i++) {
        V_out.push_back(V[v_ids[i]]);
    }
    // Fill T & others
    int cnt = 0;
    for (int i = 0; i < tetNum; i++) {
        if (t_is_removed[i]) {
            continue;
        }
        T_out.push_back(std::array<int, 4>({{map_ids[T[i][0]], map_ids[T[i][1]], map_ids[T[i][2]], map_ids[T[i][3]]}}));
        labels_out(cnt) = labels(i);
        is_surface_facet_out.push_back(is_surface_facet.at(i));
        face_on_shell_out.push_back(face_on_shell.at(i));

        cnt++;
    }

    // replace
    V = V_out;
    T = T_out;
    labels = labels_out;
    is_surface_facet = is_surface_facet_out;
    face_on_shell = face_on_shell_out;

    logger().debug("stage shell output #v = {}", V_out.size());
    logger().debug("stage shell output #t = {}", T_out.size());
}


void GetTetFromPrism(
    const std::vector<tetwild::TetVertex> &VO, 
    const std::vector<std::array<int, 4>> &TO, 
    const std::vector<std::array<int, 4>> &face_on_shell,
    const int surfaceNum,
    const prism_t &prism, 
    const std::vector<tetwild::Point_3> &VI, 
    const std::map<int, int> &map_VI2VO, 
    const std::set<int> tetsOnTargetSurface,
    // std::vector<tetwild::TetVertex> &V_temp,  // why? all vertices already exist
    std::vector<std::array<int, 4>> &T_temp, 
    std::vector<std::array<int, 4>> &is_surface_facet_temp,
    std::vector<std::array<int, 4>> &face_on_shell_temp) {

    // insert the first two tet (no further split)
    T_temp.push_back(std::array<int, 4>({{map_VI2VO.at(prism[0]), map_VI2VO.at(prism[4]), map_VI2VO.at(prism[5]), map_VI2VO.at(prism[3])}}));
    MakeTetPositive(VO, T_temp[T_temp.size()-1]);  // although these two should be OK
    T_temp.push_back(std::array<int, 4>({{map_VI2VO.at(prism[1]), map_VI2VO.at(prism[4]), map_VI2VO.at(prism[5]), map_VI2VO.at(prism[0])}}));
    MakeTetPositive(VO, T_temp[T_temp.size()-1]);
    // update their attributes
    is_surface_facet_temp.push_back(std::array<int, 4>({{-1, 1024, 1024, 1024}}));  // TODO: is this correct?
    is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1024, 1024, 1024}}));
    if (surfaceNum == 1)
        face_on_shell_temp.push_back(std::array<int, 4>({{2, 0, 0, 0}}));
    else
        face_on_shell_temp.push_back(std::array<int, 4>({{3, 0, 0, 0}}));
    face_on_shell_temp.push_back(std::array<int, 4>({{0, 0, 0, 0}}));

    // insert the third one (VO might have split one of its face to many triangles)
    Point_3 base_pt1 = VI[prism[0]];
    Point_3 base_pt2 = VI[prism[1]];
    Point_3 base_pt3 = VI[prism[2]];  // the first three vertices are the annoying triangle. 3, 4, 5 have been handled
    for (auto it=tetsOnTargetSurface.begin(); it!=tetsOnTargetSurface.end(); it++) {

        // consider the "TO_row"-th tet in TO
        int TO_row = *it;
        // "tetsOnTargetSurface" has checked whether the tet has been removed by "t_is_removed"
        for (int i=0; i<4; i++) {
            if (face_on_shell[TO_row][i] == surfaceNum) {

                Point_3 pt1 = VO[TO[TO_row][(i+1)%4]].pos;
                Point_3 pt2 = VO[TO[TO_row][(i+2)%4]].pos;
                Point_3 pt3 = VO[TO[TO_row][(i+3)%4]].pos;

                if (IsTriInsideTri(pt1, pt2, pt3, base_pt1, base_pt2, base_pt3)) {
                    // ahaha, we found a desired tet in VO. pt1, pt2, pt3 will be used as new base
                    T_temp.push_back(std::array<int, 4>({{TO[TO_row][(i+1)%4], TO[TO_row][(i+2)%4], TO[TO_row][(i+3)%4], map_VI2VO.at(prism[5])}}));

                    // update attributes
                    is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1024, 1024, 1}}));
                    face_on_shell_temp.push_back(std::array<int, 4>({{0, 0, 0, surfaceNum}}));
                    // inverted tet?
                    MakeTetPositive(VO, T_temp[T_temp.size()-1], is_surface_facet_temp[is_surface_facet_temp.size()-1], face_on_shell_temp[face_on_shell_temp.size()-1]);
                }

                // there can be at most one face of one tet that belongs to "surfaceNum"
                // because the faces of one tet cannot be colinear
                break;
            }
        }
    }
}


void GenTetMeshFromShell(
    const std::vector<tetwild::TetVertex> &VO, 
    const std::vector<std::array<int, 4>> &TO, 
    const DualShell_t &dualShell, 
    const std::vector<std::array<int, 4>> &face_on_shell,
    const std::string shellName,
    // std::vector<tetwild::TetVertex> &V_temp,  // why? all vertices already exist
    std::vector<std::array<int, 4>> &T_temp, 
    Eigen::VectorXi &labels_temp,
    std::vector<std::array<int, 4>> &is_surface_facet_temp,
    std::vector<std::array<int, 4>> &face_on_shell_temp, 
    std::vector<bool> &t_is_removed) {

    const int N = dualShell.shell_inner_bottom.size();
    int surfaceNum;
    T_temp.clear();
    is_surface_facet_temp.clear();
    face_on_shell_temp.clear();

    // Find the index of dualShell vertex in VO now
    std::map<int, int> map_VI2VO;
    if (shellName == "shell_inner_bottom") {
        MapIndex(VO, dualShell, dualShell.shell_inner_bottom, map_VI2VO);
        surfaceNum = 1;
    } else if (shellName == "shell_top_outer") {
        MapIndex(VO, dualShell, dualShell.shell_top_outer, map_VI2VO);
        surfaceNum = 4;
    }

    // Collect the row number of TO that touches "surfaceNum" surface
    std::set<int> tetsOnTargetSurface;
    for (int i=0; i<face_on_shell.size(); i++) {
        if (t_is_removed[i])
            continue;  // do not consider removed tet
        if (face_on_shell[i][0] == surfaceNum || face_on_shell[i][1] == surfaceNum || 
            face_on_shell[i][2] == surfaceNum || face_on_shell[i][3] == surfaceNum) {
                tetsOnTargetSurface.insert(i);
            }
    }

    // Generate tet mesh for each prism
    for (int i=0; i<N; i++) {

        prism_t prism;
        if (shellName == "shell_inner_bottom") {
            const prism_t &pr = dualShell.shell_inner_bottom[i];
            prism = {pr[0], pr[1], pr[2], pr[3], pr[4], pr[5]};
        } else if (shellName == "shell_top_outer") {
            const prism_t &pr = dualShell.shell_top_outer[i];
            prism = {pr[3], pr[4], pr[5], pr[0], pr[1], pr[2]};
            // Why reverse the order?
            // because i want to make the surface touching fine tets the first one
        }

        // get tet from prism & update is_surface_facet_temp + face_on_shell_temp
        GetTetFromPrism(VO, TO, face_on_shell, surfaceNum, prism, dualShell.V, map_VI2VO, tetsOnTargetSurface, T_temp, is_surface_facet_temp, face_on_shell_temp);
    }

    // we also want to enlarge "t_is_removed" accordingly
    // and of course no tet should be discarded
    std::vector<bool> falseVector(T_temp.size(), false);
    t_is_removed.insert( t_is_removed.end(), falseVector.begin(), falseVector.end() );

    // Update labels
    if (surfaceNum == 1)
        labels_temp = Eigen::VectorXi::Ones(T_temp.size(), 1);
    else if (surfaceNum == 4)
        labels_temp = Eigen::VectorXi::Ones(T_temp.size(), 1).array() * 3;
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    // now loop over all tets and label them
    // dumb way, not optimized
    for (int i=0; i<Ntet; i++) {

        Point_3 center = GetBarycenter(VO, TO, i);
        if (PointInShell(center, dualShell.shell_inner_bottom, VI))
            labels(i) = 1;
        else if (PointInShell(center, dualShell.shell_bottom_top, VI))
            labels(i) = 2;
        else if (PointInShell(center, dualShell.shell_top_outer, VI))
            labels(i) = 3;
    }

    logger().info("Tet label done");
}


void ReplaceWithPrismTet(
    const DualShell_t &dualShell, 
    std::vector<tetwild::TetVertex> &VO,
    std::vector<std::array<int, 4>> &TO, 
    Eigen::VectorXi &labels,
    std::vector<std::array<int, 4>> &is_surface_facet, 
    std::vector<std::array<int, 4>> &face_on_shell) {

    const int numTet = TO.size();

    // Remove labelled tets
    std::vector<bool> t_is_removed(TO.size(), false);
    int numOldTet = 0;
    for (int i=0; i<numTet; i++) {
        if (labels(i) != 0) 
            t_is_removed[i] = true;
        else
            numOldTet++;
    }
    logger().debug("#tets after removing labeled shell tets = {}", numOldTet);

    // Generate new tets
    std::vector<std::array<int, 4>> T_temp;
    std::vector<std::array<int, 4>> is_surface_facet_temp;
    std::vector<std::array<int, 4>> face_on_shell_temp;
    Eigen::VectorXi labels_temp;
    /////// FOR shell_inner_bottom
    GenTetMeshFromShell(VO, TO, dualShell, face_on_shell, "shell_inner_bottom", T_temp, labels_temp, is_surface_facet_temp, face_on_shell_temp, t_is_removed);
    // concatenate new tet with old
    TO.insert( TO.end(), T_temp.begin(), T_temp.end() );
    is_surface_facet.insert( is_surface_facet.end(), is_surface_facet_temp.begin(), is_surface_facet_temp.end() );
    face_on_shell.insert( face_on_shell.end(), face_on_shell_temp.begin(), face_on_shell_temp.end() );
      Eigen::VectorXi temp = labels;
      labels.resize(temp.rows() + labels_temp.rows());
      labels << temp, labels_temp;
    logger().debug("SHELL_INNER_BOTTOM done");
    ////// FOR shell_top_outer
    GenTetMeshFromShell(VO, TO, dualShell, face_on_shell, "shell_top_outer", T_temp, labels_temp, is_surface_facet_temp, face_on_shell_temp, t_is_removed);
    // concatenate new tet with old
    TO.insert( TO.end(), T_temp.begin(), T_temp.end() );
    is_surface_facet.insert( is_surface_facet.end(), is_surface_facet_temp.begin(), is_surface_facet_temp.end() );
    face_on_shell.insert( face_on_shell.end(), face_on_shell_temp.begin(), face_on_shell_temp.end() );
      temp = labels;
      labels.resize(temp.rows() + labels_temp.rows());
      labels << temp, labels_temp;
    logger().debug("SHELL_TOP_OUTER done");

    // remove inplace
    CleanTetMesh(t_is_removed, VO, TO, labels, is_surface_facet, face_on_shell);
    std::cout << labels << std::endl;
    logger().info("Replace with prism tet done");
}

}  // namespace tetshell
