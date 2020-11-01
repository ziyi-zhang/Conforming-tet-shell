#include <shell/common.hpp>
#include <shell/Shell.h>
#include <shell/Utils.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>

#include <Eigen/Dense>


namespace tetshell {

using namespace std;
using tetwild::Point_3;
using tetwild::Segment_3;
using tetwild::logger;

namespace {

void ConstructShellFromTri(const RowMatX3i &F, int M, int a, int b, shell_t &shell) {

    const int N = F.rows();
    if (a<0 || b<0 || a>3 || b>3 || a==b)
        tetwild::log_and_throw("ConstructShellFromTri invalid input");

    shell.clear();
    for (int i=0; i<N; i++) {

        prism_t prism{F(i, 0) + a*M, F(i, 1) + a*M, F(i, 2) + a*M, 
                      F(i, 0) + b*M, F(i, 1) + b*M, F(i, 2) + b*M};
        /*
        // we no longer store coordinates here. Only store the index in VI
        prism[0] = V_bottom[F(i, 0)];
        prism[1] = V_bottom[F(i, 1)];
        prism[2] = V_bottom[F(i, 2)];
        prism[3] = V_top   [F(i, 0)];
        prism[4] = V_top   [F(i, 1)];
        prism[5] = V_top   [F(i, 2)];
        */
        shell.push_back(prism);
    }
}


bool IsDegeneratedTet(const Point_3 &pt1, const Point_3 &pt2, const Point_3 &pt3, const Point_3 &pt4) {

    return (pt1 == pt2) || (pt1 == pt3) || (pt1 == pt4) || (pt2 == pt3) || (pt2 == pt4) || (pt3 == pt4);
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


/// Is point pt1 inside triangle (base_pt0, base_pt1, base_pt2)
bool IsPointInsideTri(
    const Point_3 &pt1, 
    const Point_3 &base_pt0, 
    const Point_3 &base_pt1, 
    const Point_3 &base_pt2) {

    tetwild::Triangle_3 tri(base_pt0, base_pt1, base_pt2);
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
    for (int i = 0; i < T.size(); i++) {
        if (t_is_removed[i]) {
            continue;
        }
        T_out.push_back(std::array<int, 4>({{map_ids[T[i][0]], map_ids[T[i][1]], map_ids[T[i][2]], map_ids[T[i][3]]}}));
        labels_out(cnt) = labels(i);
        is_surface_facet_out.push_back(is_surface_facet.at(i));
        face_on_shell_out.push_back(face_on_shell.at(i));

        cnt++;
    }

    // back substitute
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
    const std::vector<std::array<int, 4>> &face_on_shell,  // The old "face_on_shell", recording which input surface each face of a tet is on
    const int surfaceIdx,  // either "SURFACE_INNER" or "SURFACE_OUTER", the "bad" subdivided face
    const prism_t &prism,  // the vertex index in "prism" is w.r.t. VI, so we need "map_VI2VO"
    const std::vector<tetwild::Point_3> &VI, 
    const std::map<int, int> &map_VI2VO, 
    const std::vector<bool> &t_is_removed,
    // std::vector<tetwild::TetVertex> &V_temp,  // why? all vertices already exist
    std::vector<std::array<int, 4>> &T_temp, 
    std::vector<std::array<int, 4>> &is_surface_facet_temp,
    std::vector<std::array<int, 4>> &face_on_shell_temp) {

    // This is the prism
    ///   3 ------- 5
    ///   | \__4__/ |
    ///   |    |    |
    ///   0----|--- 2
    ///     \__1__/
    // We split it into three tet {0, 4, 5, 3} {1, 4, 5, 0} and {0, 1, 2, 5}
    // It is guaranteed that triangle (345) is intact and only face (012) may be subdivided

    // Define the bottom triangle (012)
    Point_3 base_pt0 = VI[prism[0]];
    Point_3 base_pt1 = VI[prism[1]];
    Point_3 base_pt2 = VI[prism[2]];  // the annoying triangle
    Segment_3 edge_01(base_pt0, base_pt1);
    // Collect vertices in VO that are on edge (01)
    // sort "vertsOnTargetEdge" in order from "base_pt0" to "base_pt1"
    std::map<tetwild::CGAL_FT, int> vertsOnTargetEdge;
    // Collect the row number of TO that touches "surfaceIdx" surface
    std::unordered_set<int> tetsOnTargetSurface;  // TO index
    // fill the two unordered_sets
    for (int i=0; i<face_on_shell.size(); i++) {
        if (t_is_removed[i])
            continue;  // do not consider removed tet
        if (face_on_shell[i][0] == surfaceIdx || face_on_shell[i][1] == surfaceIdx || 
            face_on_shell[i][2] == surfaceIdx || face_on_shell[i][3] == surfaceIdx) {
                // This tet is in
                tetsOnTargetSurface.insert(i);
                // further check edge (01)
                for (int j=0; j<4; j++) {
                    const Point_3 &pt = VO[TO[i][j]].pos;
                    if (edge_01.has_on(pt)) {
                        // "pt" is on edge_01, insert to set
                        Segment_3 edge_temp(base_pt0, pt);
                        vertsOnTargetEdge.insert(std::make_pair(edge_temp.squared_length(), TO[i][j]));  // will be sorted
                    }
                }
            }   
    }

    // Stats of degenerated tet
    int cnt_singularity_type1 = 0;
    int cnt_singularity_type2 = 0;
    int cnt_singularity_type3 = 0;

    ///////////////
    //   FIRST   //
    ///////////////
    // insert the first tet {0, 4, 5, 3} (no further split)
    if (!IsDegeneratedTet(VI[prism[0]], VI[prism[4]], VI[prism[5]], VI[prism[3]])) {  // no singularity

        T_temp.push_back(std::array<int, 4>({{map_VI2VO.at(prism[0]), map_VI2VO.at(prism[4]), map_VI2VO.at(prism[5]), map_VI2VO.at(prism[3])}}));
        // update this new tet's attribute
        if (IsTetPositive(VI[prism[4]], VI[prism[5]], VI[prism[3]], VI[prism[0]]))
            is_surface_facet_temp.push_back(std::array<int, 4>({{1, 1024, 1024, 1024}}));
        else
            is_surface_facet_temp.push_back(std::array<int, 4>({{-1, 1024, 1024, 1024}}));
        if (surfaceIdx == SURFACE_INNER)
            face_on_shell_temp.push_back(std::array<int, 4>({{SURFACE_BOTTOM, NOT_SUR, NOT_SUR, NOT_SUR}}));
        else if (surfaceIdx == SURFACE_OUTER)
            face_on_shell_temp.push_back(std::array<int, 4>({{SURFACE_TOP, NOT_SUR, NOT_SUR, NOT_SUR}}));
        else
            tetwild::log_and_throw("GetTetFromPrism: surfaceIdx invalid");
        // positive tet
        if (MakeTetPositive(VO, T_temp.back(), is_surface_facet_temp.back(), face_on_shell_temp.back())) {
            is_surface_facet_temp.back()[2] *= -1;
            /*
            // DEBUG PURPOSE
            std::array<double, 3> pt_1 = {CGAL::to_double(VI[prism[0]][0]), CGAL::to_double(VI[prism[0]][1]), CGAL::to_double(VI[prism[0]][2])};
            std::array<double, 3> pt_2 = {CGAL::to_double(VI[prism[1]][0]), CGAL::to_double(VI[prism[1]][1]), CGAL::to_double(VI[prism[1]][2])};
            std::array<double, 3> pt_3 = {CGAL::to_double(VI[prism[2]][0]), CGAL::to_double(VI[prism[2]][1]), CGAL::to_double(VI[prism[2]][2])};
            std::array<double, 3> pt_4 = {CGAL::to_double(VI[prism[3]][0]), CGAL::to_double(VI[prism[3]][1]), CGAL::to_double(VI[prism[3]][2])};
            std::array<double, 3> pt_5 = {CGAL::to_double(VI[prism[4]][0]), CGAL::to_double(VI[prism[4]][1]), CGAL::to_double(VI[prism[4]][2])};
            std::array<double, 3> pt_6 = {CGAL::to_double(VI[prism[5]][0]), CGAL::to_double(VI[prism[5]][1]), CGAL::to_double(VI[prism[5]][2])};
            std::cerr << "first tet being inverted" << std::endl;
            */
        }
    } else {
        cnt_singularity_type1++;
    }

    ////////////////
    //   SECOND   //
    ////////////////
    if (!IsDegeneratedTet(VI[prism[1]], VI[prism[4]], VI[prism[5]], VI[prism[0]])) {  // no singularity

        // check the first and last of "vertsOnTargetEdge" is pt0 and pt1
        if (base_pt0 != VO[vertsOnTargetEdge.begin()->second].pos || base_pt1 != VO[vertsOnTargetEdge.rbegin()->second].pos) {
            tetwild::log_and_throw("GetTetFromPrism: vertsOnTargetEdge construction error.");
        }
        // insert the second tet {1, 4, 5, 0} (some split due to subdivision on edge (01))
        int ptIdx_1, ptIdx_2;
        const int ptIdx_3 = map_VI2VO.at(prism[5]);  // the two vertices from the prism
        const int ptIdx_4 = map_VI2VO.at(prism[4]);  // the two vertices from the prism
        for (auto it=vertsOnTargetEdge.begin(); it!=vertsOnTargetEdge.end();) {

            ptIdx_1 = it->second;
            it++;
            if (it == vertsOnTargetEdge.end()) break;
            ptIdx_2 = it->second;

            T_temp.push_back(std::array<int, 4>({{ptIdx_1, ptIdx_2, ptIdx_3, ptIdx_4}}));
            // update this new tet's attribute
            is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1024, 1024, 1024}}));
            face_on_shell_temp.push_back(std::array<int, 4>({{NOT_SUR, NOT_SUR, NOT_SUR, NOT_SUR}}));
            // positive tet
            MakeTetPositive(VO, T_temp.back(), is_surface_facet_temp.back(), face_on_shell_temp.back());
        }
    } else {
        cnt_singularity_type2++;
    }

    ///////////////
    //   THIRD   //
    ///////////////
    // insert the third one (VO might have split face (012) to many triangles)
    if (!IsDegeneratedTet(VI[prism[0]], VI[prism[1]], VI[prism[2]], VI[prism[5]])) {  // no singularity

        for (auto it=tetsOnTargetSurface.begin(); it!=tetsOnTargetSurface.end(); it++) {

            // consider the "TO_row"-th tet in TO
            int TO_row = *it;
            int ptIdx1, ptIdx2, ptIdx3;
            const int ptIdx4 = map_VI2VO.at(prism[5]);  // the vertex from the prism
            // "tetsOnTargetSurface" has checked whether the tet has been removed by "t_is_removed"
            for (int i=0; i<4; i++) {
                if (face_on_shell[TO_row][i] == surfaceIdx) {

                    ptIdx1 = TO[TO_row][(i+1)%4];
                    ptIdx2 = TO[TO_row][(i+2)%4];
                    ptIdx3 = TO[TO_row][(i+3)%4];
                    Point_3 pt1 = VO[ptIdx1].pos;
                    Point_3 pt2 = VO[ptIdx2].pos;
                    Point_3 pt3 = VO[ptIdx3].pos;

                    if (IsTriInsideTri(pt1, pt2, pt3, base_pt0, base_pt1, base_pt2)) {
                        // ahaha, we found a desired tet in VO. pt1, pt2, pt3 will be used as new base
                        T_temp.push_back(std::array<int, 4>({{ptIdx1, ptIdx2, ptIdx3, ptIdx4}}));
                        // update attributes
                        if (IsTetPositive(pt1, pt2, pt3, VI[prism[5]]))
                            is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1024, 1024, 1}}));
                        else
                            is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1024, 1024, -1}}));
                        face_on_shell_temp.push_back(std::array<int, 4>({{0, 0, 0, surfaceIdx}}));
                        // positive tet
                        if (MakeTetPositive(VO, T_temp.back(), is_surface_facet_temp.back(), face_on_shell_temp.back())) {
                            is_surface_facet_temp.back()[3] *= -1;
                        }
                    }

                    // there can be at most one face of any tet that is on "surfaceIdx"
                    // because the faces of one tet cannot be colinear
                    break;
                }
            }
        }
    } else {
        cnt_singularity_type3++;
    }

    
}


void GenTetMeshFromShell(
    const std::vector<tetwild::TetVertex> &VO, 
    const std::vector<std::array<int, 4>> &TO, 
    const DualShell_t &dualShell, 
    const std::vector<std::array<int, 4>> &face_on_shell,
    const int shellName,  // either SHELL_INNER_BOTTOM or SHELL_TOP_OUTER
    // std::vector<tetwild::TetVertex> &V_temp,  // why? all vertices already exist
    std::vector<std::array<int, 4>> &T_temp, 
    Eigen::VectorXi &labels_temp,
    std::vector<std::array<int, 4>> &is_surface_facet_temp,
    std::vector<std::array<int, 4>> &face_on_shell_temp, 
    std::vector<bool> &t_is_removed) {

    const int N = dualShell.shell_inner_bottom.size();
    int surfaceIdx;
    T_temp.clear();
    is_surface_facet_temp.clear();
    face_on_shell_temp.clear();

    // Find the index of dualShell vertex in VO now
    std::map<int, int> map_VI2VO;
    if (shellName == SHELL_INNER_BOTTOM) {
        MapIndex(VO, dualShell, dualShell.shell_inner_bottom, map_VI2VO);
        surfaceIdx = SURFACE_INNER;  // the one of the two surfaces that has been subdivided 
    } else if (shellName == SHELL_TOP_OUTER) {
        MapIndex(VO, dualShell, dualShell.shell_top_outer, map_VI2VO);
        surfaceIdx = SURFACE_OUTER;  // the one of the two surfaces that has been subdivided 
    }

    // Generate tet mesh for each prism
    for (int i=0; i<N; i++) {

        prism_t prism;
        if (shellName == SHELL_INNER_BOTTOM) {
            const prism_t &pr = dualShell.shell_inner_bottom[i];
            prism = {pr[0], pr[1], pr[2], pr[3], pr[4], pr[5]};
        } else if (shellName == SHELL_TOP_OUTER) {
            const prism_t &pr = dualShell.shell_top_outer[i];
            prism = {pr[3], pr[4], pr[5], pr[0], pr[1], pr[2]};
            // Why reverse the order?
            // to make the surface touching fine tets always the first one
        }

        // get tet from prism & update is_surface_facet_temp + face_on_shell_temp
        GetTetFromPrism(VO, TO, face_on_shell, surfaceIdx, prism, dualShell.V, map_VI2VO, t_is_removed,  // const input
                        T_temp, is_surface_facet_temp, face_on_shell_temp);  // output
    }

    // we also want to enlarge "t_is_removed" accordingly
    // and of course no tet should be discarded
    std::vector<bool> falseVector(T_temp.size(), false);
    t_is_removed.insert( t_is_removed.end(), falseVector.begin(), falseVector.end() );

    // Update labels
    // All the tets generated here must be in the same region
    if (surfaceIdx == SURFACE_INNER)
        labels_temp = Eigen::VectorXi::Ones(T_temp.size(), 1).array() * SHELL_INNER_BOTTOM;
    else if (surfaceIdx == SURFACE_OUTER)
        labels_temp = Eigen::VectorXi::Ones(T_temp.size(), 1).array() * SHELL_TOP_OUTER;
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////

void GenDualShell(const std::vector<Point_3> &VI, const Eigen::MatrixXi &FI, DualShell_t &dualShell) {

    const int M = VI.size() / 4;
    const int Ntri = FI.rows() / 4;

    // assert size - check this even in release
    if ((VI.size() % 4 != 0) ||
        (FI.rows() % 4 != 0)) {
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
        tetwild::log_and_throw("GenDualShell: Wrong input index matrix");
    }

    // split input
    dualShell.F = F1;
    dualShell.V = VI;
    /*
    dualShell.V_inner = {VI.begin(), VI.begin()+M};
    dualShell.V_bottom = {VI.begin()+M, VI.begin()+2*M};
    dualShell.V_top = {VI.begin()+2*M, VI.begin()+3*M};
    dualShell.V_outer = {VI.begin()+3*M, VI.begin()+4*M};
    */

    // Retrieve prisms from dual-shell
    ConstructShellFromTri(F1, M, 0, 1, dualShell.shell_inner_bottom);
    ConstructShellFromTri(F1, M, 1, 2, dualShell.shell_bottom_top);
    ConstructShellFromTri(F1, M, 2, 3, dualShell.shell_top_outer);
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

    //// FOR shell_inner_bottom
    GenTetMeshFromShell(VO, TO, dualShell, face_on_shell, SHELL_INNER_BOTTOM, T_temp, labels_temp, is_surface_facet_temp, face_on_shell_temp, t_is_removed);
    // concatenate new tet with old
    TO.insert( TO.end(), T_temp.begin(), T_temp.end() );
    is_surface_facet.insert( is_surface_facet.end(), is_surface_facet_temp.begin(), is_surface_facet_temp.end() );
    face_on_shell.insert( face_on_shell.end(), face_on_shell_temp.begin(), face_on_shell_temp.end() );
      Eigen::VectorXi temp = labels;
      labels.resize(temp.rows() + labels_temp.rows());
      labels << temp, labels_temp;
    logger().debug("GenTetMeshFromShell: SHELL_INNER_BOTTOM done");

    /// FOR shell_top_outer
    GenTetMeshFromShell(VO, TO, dualShell, face_on_shell, SHELL_TOP_OUTER, T_temp, labels_temp, is_surface_facet_temp, face_on_shell_temp, t_is_removed);
    // concatenate new tet with old
    TO.insert( TO.end(), T_temp.begin(), T_temp.end() );
    is_surface_facet.insert( is_surface_facet.end(), is_surface_facet_temp.begin(), is_surface_facet_temp.end() );
    face_on_shell.insert( face_on_shell.end(), face_on_shell_temp.begin(), face_on_shell_temp.end() );
      temp = labels;
      labels.resize(temp.rows() + labels_temp.rows());
      labels << temp, labels_temp;
    logger().debug("GenTetMeshFromShell: SHELL_TOP_OUTER done");

    // remove "t_is_removed" inplace
    CleanTetMesh(t_is_removed, VO, TO, labels, is_surface_facet, face_on_shell);
    logger().info("Replace with prism tet done");
}

}  // namespace tetshell
