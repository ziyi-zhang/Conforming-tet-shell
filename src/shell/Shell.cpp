#include <shell/common.hpp>
#include <shell/Shell.h>
#include <shell/Utils.h>
#include <tetwild/Args.h>
#include <tetwild/State.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>
#include <tetwild/TetmeshElements.h>

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
                std::cerr << "VI_index = " << VI_index << std::endl;
                std::cerr << "V = " << CGAL::to_double(VI[VI_index][0]) << " " << CGAL::to_double(VI[VI_index][1]) << " " <<  CGAL::to_double(VI[VI_index][2]) << std::endl;
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

    logger().debug("stage shell output #VO = {}", V_out.size());
    logger().debug("stage shell output #TO = {}", T_out.size());
}


void GetTetFromPrism(
    const std::vector<tetwild::TetVertex> &VO, 
    const std::vector<std::array<int, 4>> &TO, 
    const std::unordered_set<int> &tetsOnTargetSurface,  // pre-calculated tetsOnTargetSurface
    const std::unordered_set<int> &vertsOnTargetEdge_candidate,  // pre-calculated vertsOnTargetEdge_candidate
    const std::vector<std::array<int, 4>> &face_on_shell,  // The old "face_on_shell", recording which input surface each face of a tet is on
    const int surfaceIdx,  // either "SURFACE_INNER" or "SURFACE_OUTER", the "bad" subdivided face
    const prism_t &prism,  // the vertex index in "prism" is w.r.t. VI, so we need "map_VI2VO"
    const std::vector<tetwild::Point_3> &VI, 
    const std::map<int, int> &map_VI2VO, 
    const std::vector<bool> &t_is_removed,
    // std::vector<tetwild::TetVertex> &V_temp,  // why? all vertices already exist
    std::vector<std::array<int, 4>> &T_temp, 
    std::vector<std::array<int, 4>> &is_surface_facet_temp,
    std::vector<std::array<int, 4>> &face_on_shell_temp,
    std::vector<int> &labels_temp,  // DEBUG PURPOSE, should be filled with the same number
    std::array<int, 3> &singularCnt) {

    // This is the prism (not necessarily a pentahedron)
    ///   3 ------- 5
    ///   | \__4__/ |
    ///   |    |    |
    ///   0----|--- 2
    ///     \__1__/
    // We split it into three tet 
    //      if idx(1)>idx(2): {0, 3, 4, 5}, {0, 5, 4, 2} and {0, 1, 2, 4}
    //      if idx(1)<idx(2): {0, 3, 4, 5}, {0, 5, 4, 1} and {0, 1, 2, 5}
    // It is guaranteed that triangle (345) is intact and only face (012) may be subdivided

    // singualrity
    bool top_collapse;  // whether the top has collapsed due to singularity
    bool middle_collapse;
    bool bottom_collapse;  // whether bottom has collapsed due to singularity
    // Define the bottom triangle (012)
    Point_3 base_pt0 = VI[prism[0]];
    Point_3 base_pt1 = VI[prism[1]];
    Point_3 base_pt2 = VI[prism[2]];  // the annoying triangle
    // check
    if (!(prism[0] < prism[1] && prism[0] < prism[2])) {
        tetwild::log_and_throw("GetTetFromPrism: Invalid prism. prism[0] not the smallest index.");
    }
    // consistent split method
    // NOTE: comparison is based on VI-index, not VO-index
    bool tetSplitA = (prism[1] > prism[2]);
    auto tetSplit = tetSplitA ? TETRA_SPLIT_A : TETRA_SPLIT_B;
    Point_3 &base_ptx = tetSplitA ? base_pt2 : base_pt1;  // 'x' is either 1 or 2
    Segment_3 edge_0x(base_pt0, base_ptx);
    /*
    // NOTE: This code is used for debug. LOW efficient
    // Collect vertices in VO that are on edge (0x)
    // sort "vertsOnTargetEdge" in order from "base_pt0" to "base_ptx"
    std::map<tetwild::CGAL_FT, int> vertsOnTargetEdge;
    // Collect the row number of TO that touches "surfaceIdx" surface
    std::unordered_set<int> tetsOnTargetSurface;  // TO index
    // fill "vertsOnTargetEdge" and "tetsOnTargetSurface"
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
                    if (edge_0x.has_on(pt)) {
                        // "pt" is on edge_0x, insert to set
                        Segment_3 edge_temp(base_pt0, pt);
                        vertsOnTargetEdge.insert(std::make_pair(edge_temp.squared_length(), TO[i][j]));  // will be sorted
                    }
                }
            }   
    }
    */
    // Collect vertices in VO that are on edge (0x)
    // sort "vertsOnTargetEdge" in order from "base_pt0" to "base_ptx"
    std::map<tetwild::CGAL_FT, int> vertsOnTargetEdge;
    for (auto it=vertsOnTargetEdge_candidate.begin(); it!=vertsOnTargetEdge_candidate.end(); it++) {
        const Point_3 &pt = VO[*it].pos;
        if (edge_0x.has_on(pt)) {
            // "pt" is on edge_0x, insert to set
            Segment_3 edge_temp(base_pt0, pt);
            vertsOnTargetEdge.insert(std::make_pair(edge_temp.squared_length(), *it));  // will be sorted
        }
    }

    /////////////////
    // SINGULARITY //
    /////////////////
    // TYPE - 1
    top_collapse = IsDegeneratedTet(VI[prism[0]], VI[prism[3]], VI[prism[4]], VI[prism[5]]);
    // TYPE - 2
    //     (prism[1] > prism[2])?        True : False
    int local_idx1 = tetSplit[1][0];  //  0       0
    int local_idx2 = tetSplit[1][1];  //  5       5
    int local_idx3 = tetSplit[1][2];  //  4       4
    int local_idx4 = tetSplit[1][3];  //  2       1  (== x)
    middle_collapse = IsDegeneratedTet(VI[prism[local_idx1]], VI[prism[local_idx2]], VI[prism[local_idx3]], VI[prism[local_idx4]]);
    // TYPE - 3
    const Point_3 &pt_tip = VI[prism[tetSplit[2][3]]];  // 4 or 5
    bottom_collapse = IsDegeneratedTet(VI[prism[0]], VI[prism[1]], VI[prism[2]], pt_tip);
    // assert
    if ((top_collapse && middle_collapse) || (top_collapse && bottom_collapse) || (middle_collapse && bottom_collapse)) {
        tetwild::log_and_throw("GetTetFromPrism: adjacent singularities found.");
    }

    ////////////////
    //   TYPE 1   //
    ////////////////
    // insert the first tet {0, 3, 4, 5} (no further split)
    if (!top_collapse) {  // no singularity

        T_temp.push_back(std::array<int, 4>({{map_VI2VO.at(prism[0]), map_VI2VO.at(prism[3]), map_VI2VO.at(prism[4]), map_VI2VO.at(prism[5])}}));
        // update this new tet's attribute
        labels_temp.push_back(7);  // DEBUG PURPOSE
        is_surface_facet_temp.push_back(std::array<int, 4>({{1, 1024, 1024, 1024}}));  // force it to be 1
        if (surfaceIdx == SURFACE_INNER)
            face_on_shell_temp.push_back(std::array<int, 4>({{SURFACE_BOTTOM, NOT_SUR, NOT_SUR, NOT_SUR}}));
        else if (surfaceIdx == SURFACE_OUTER)
            face_on_shell_temp.push_back(std::array<int, 4>({{SURFACE_TOP, NOT_SUR, NOT_SUR, NOT_SUR}}));
        else
            tetwild::log_and_throw("GetTetFromPrism: surfaceIdx invalid");
        // positive tet
        if (MakeTetPositive(VO, T_temp.back(), is_surface_facet_temp.back(), face_on_shell_temp.back())) {
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
        // logger().warn("{}, {}, {}, {}, {}, {}", prism[0], prism[1], prism[2], prism[3], prism[4], prism[5]);
        singularCnt[0]++;
    }

    ////////////////
    //   TYPE 2   //
    ////////////////
    if (!middle_collapse) {  // no singularity

        // check the first and last of "vertsOnTargetEdge" is pt0 and ptx
        if (base_pt0 != VO[vertsOnTargetEdge.begin()->second].pos || base_ptx != VO[vertsOnTargetEdge.rbegin()->second].pos) {
            
            printf("vertsOnTargetEdge.size = %lu\n", vertsOnTargetEdge.size());
            PrintPoints(base_pt0, VO[vertsOnTargetEdge.begin()->second].pos, base_ptx, VO[vertsOnTargetEdge.rbegin()->second].pos);
            tetwild::log_and_throw("GetTetFromPrism: vertsOnTargetEdge construction error.");
        }
        // insert the second tet {0, 5, 4, 2}/{0, 5, 4, 1} (some split due to subdivision on edge (0x))
        int ptIdx1, ptIdx4;
        const int ptIdx2 = map_VI2VO.at(prism[5]);  // the two vertices from the prism
        const int ptIdx3 = map_VI2VO.at(prism[4]);  // the two vertices from the prism
        for (auto it=vertsOnTargetEdge.begin(); it!=vertsOnTargetEdge.end();) {

            bool is_tet_with_idx0 = (it == vertsOnTargetEdge.begin());
            ptIdx1 = it->second;
            it++;
            if (it == vertsOnTargetEdge.end()) break;
            ptIdx4 = it->second;

            T_temp.push_back(std::array<int, 4>({{ptIdx1, ptIdx2, ptIdx3, ptIdx4}}));
            // update this new tet's attribute
            labels_temp.push_back(8);  // DEBUG PURPOSE
            /// NOTE: the face label depends on whether some collapse has occured
            if (top_collapse && is_tet_with_idx0) {  // this correction only works for the tet with face 045 (not for the subdivided ones)
                // top has collapsed
                int surface_label;
                if (surfaceIdx == SURFACE_INNER) {
                    surface_label = SURFACE_BOTTOM;
                } else if (surfaceIdx == SURFACE_OUTER) {
                    surface_label = SURFACE_TOP;
                } else {
                    tetwild::log_and_throw("surfaceIdx error.");
                }
                is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1024, 1024, 1}}));
                face_on_shell_temp.push_back(std::array<int, 4>({{NOT_SUR, NOT_SUR, NOT_SUR, surface_label}}));
            } else if (bottom_collapse) {  // this correction works for all subdivided tets
                // bottom has collapsed
                if (tetSplitA) {
                    is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1, 1024, 1024}}));
                    face_on_shell_temp.push_back(std::array<int, 4>({{NOT_SUR, surfaceIdx, NOT_SUR, NOT_SUR}}));
                } else {
                    is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1024, 1, 1024}}));
                    face_on_shell_temp.push_back(std::array<int, 4>({{NOT_SUR, NOT_SUR, surfaceIdx, NOT_SUR}}));
                }
            } else {
                // ordinary
                is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1024, 1024, 1024}}));
                face_on_shell_temp.push_back(std::array<int, 4>({{NOT_SUR, NOT_SUR, NOT_SUR, NOT_SUR}}));
            }
            // positive tet
            if (MakeTetPositive(VO, T_temp.back())) {
            }
        }
    } else {
        singularCnt[1]++;
    }

    ////////////////
    //   TYPE 3   //
    ////////////////
    // insert the third one {0, 1, 2, 4}/{0, 1, 2, 5} (VO might have split face (012) to many triangles)
    int cnt_inserted = 0;
    int ptIdx1, ptIdx2, ptIdx3;
    if (!bottom_collapse) {  // no singularity

        for (auto it=tetsOnTargetSurface.begin(); it!=tetsOnTargetSurface.end(); it++) {

            // consider the "TO_row"-th tet in TO
            int TO_row = *it;
            const int ptIdx4 = map_VI2VO.at(prism[tetSplit[2][3]]);  // the vertex from the prism
            // "tetsOnTargetSurface" has checked whether the tet has been removed by "t_is_removed"
            for (int i=0; i<4; i++) {
                if (face_on_shell[TO_row][i] == surfaceIdx) {

                    ptIdx1 = TO[TO_row][(i+1)%4];
                    ptIdx2 = TO[TO_row][(i+2)%4];
                    ptIdx3 = TO[TO_row][(i+3)%4];
                    Point_3 pt1 = VO[ptIdx1].pos;
                    Point_3 pt2 = VO[ptIdx2].pos;
                    Point_3 pt3 = VO[ptIdx3].pos;

                    /*
                    // DEBUG PURPOSE
                    if (ptIdx1 == 580 || ptIdx1 == 587 || ptIdx1 == 541) {
                        int pt1_x = CGAL::to_double(pt1[0]);
                        int pt1_y = CGAL::to_double(pt1[1]);
                        int pt1_z = CGAL::to_double(pt1[2]);

                        int pt2_x = CGAL::to_double(pt2[0]);
                        int pt2_y = CGAL::to_double(pt2[1]);
                        int pt2_z = CGAL::to_double(pt2[2]);

                        int pt3_x = CGAL::to_double(pt3[0]);
                        int pt3_y = CGAL::to_double(pt3[1]);
                        int pt3_z = CGAL::to_double(pt3[2]);

                        bool pt1_good = IsPointInsideTri(pt1, base_pt0, base_pt1, base_pt2);
                        bool pt2_good = IsPointInsideTri(pt2, base_pt0, base_pt1, base_pt2);
                        bool pt3_good = IsPointInsideTri(pt3, base_pt0, base_pt1, base_pt2);

                        logger().info("{} {} {}: {} {} {}", ptIdx1, ptIdx2, ptIdx3, pt1_good, pt2_good, pt3_good);
                    }
                    */

                    if (IsTriInsideTri(pt1, pt2, pt3, base_pt0, base_pt1, base_pt2)) {
                        // aha, we found a desired tet in VO. pt1, pt2, pt3 will be used as new base
                        T_temp.push_back(std::array<int, 4>({{ptIdx1, ptIdx2, ptIdx3, ptIdx4}}));
                        // update attributes
                        cnt_inserted++;
                        labels_temp.push_back(9);  // DEBUG PURPOSE
                        is_surface_facet_temp.push_back(std::array<int, 4>({{1024, 1024, 1024, 1}}));
                        face_on_shell_temp.push_back(std::array<int, 4>({{NOT_SUR, NOT_SUR, NOT_SUR, surfaceIdx}}));
                        // positive tet
                        if (MakeTetPositive(VO, T_temp.back(), is_surface_facet_temp.back(), face_on_shell_temp.back())) {

                        }
                    }

                    // there can be more than one face of a tet that is on "surfaceIdx"
                    //// DO NOT BREAK HERE
                    // break;
                }
            }
        }
        if (cnt_inserted == 0) {
            logger().error("prism0 = {} prism1 = {} prism2 = {}", prism[0], prism[1], prism[2]);
            logger().error("VO_1 = {}, VO_2 = {}, VO_3 = {}", map_VI2VO.at(prism[0]), map_VI2VO.at(prism[1]), map_VI2VO.at(prism[2]));
            double x1 = CGAL::to_double(VO[ map_VI2VO.at(prism[0]) ].pos[0]);
            double y1 = CGAL::to_double(VO[ map_VI2VO.at(prism[0]) ].pos[1]);
            double z1 = CGAL::to_double(VO[ map_VI2VO.at(prism[0]) ].pos[2]);
            double x2 = CGAL::to_double(VO[ map_VI2VO.at(prism[1]) ].pos[0]);
            double y2 = CGAL::to_double(VO[ map_VI2VO.at(prism[1]) ].pos[1]);
            double z2 = CGAL::to_double(VO[ map_VI2VO.at(prism[1]) ].pos[2]);
            double x3 = CGAL::to_double(VO[ map_VI2VO.at(prism[2]) ].pos[0]);
            double y3 = CGAL::to_double(VO[ map_VI2VO.at(prism[2]) ].pos[1]);
            double z3 = CGAL::to_double(VO[ map_VI2VO.at(prism[2]) ].pos[2]);
            logger().error("x1={} y1={} z1={} x2={} y2={} z2={} x3={} y3={} z3={}", x1, y1, z1, x2, y2, z2, x3, y3, z3);
            tetwild::log_and_throw("GetTetFromPrism: Type 3 tetrahedron not founded.");
        }
    } else {
        singularCnt[2]++;
    }
}


void GenTetMeshFromShell(
    const tetwild::Args &args, 
    const std::vector<tetwild::TetVertex> &VO, 
    const std::vector<std::array<int, 4>> &TO, 
    const DualShell_t &dualShell, 
    const std::vector<std::array<int, 4>> &face_on_shell,
    const std::vector<std::array<int, 4>> &faceIdx_on_shell,
    const int shellName,  // either SHELL_INNER_BOTTOM or SHELL_TOP_OUTER
    // std::vector<tetwild::TetVertex> &V_temp,  // why? all vertices already exist
    std::vector<std::array<int, 4>> &T_temp, 
    Eigen::VectorXi &labels_temp,  // DEBUG purpose, should be filled with the same number
    std::vector<std::array<int, 4>> &is_surface_facet_temp,
    std::vector<std::array<int, 4>> &face_on_shell_temp, 
    std::vector<bool> &t_is_removed) {

    const int N = dualShell.shell_inner_bottom.size();
    int surfaceIdx;
    std::array<int, 3> singularCnt = {0, 0, 0};
    T_temp.clear();
    is_surface_facet_temp.clear();
    face_on_shell_temp.clear();
    std::vector<int> labels_temp_vec;  // vector version of labels_temp

    // Find the index of dualShell vertex in VO now
    std::map<int, int> map_VI2VO;
    if (shellName == SHELL_INNER_BOTTOM) {
        MapIndex(VO, dualShell, dualShell.shell_inner_bottom, map_VI2VO);
        surfaceIdx = SURFACE_INNER;  // the one of the two surfaces that has been subdivided 
    } else if (shellName == SHELL_TOP_OUTER) {
        MapIndex(VO, dualShell, dualShell.shell_top_outer, map_VI2VO);
        surfaceIdx = SURFACE_OUTER;  // the one of the two surfaces that has been subdivided 
    }
    logger().info("MapIndex calculated");

    // NOTE: for the explanation of this part, see the commented area in "GetTetFromPrism"
    // Calculate two data structures in advance
    std::vector<std::unordered_set<int> > vertsOnTargetEdge_candidate_array(N);
    std::vector<std::unordered_set<int> > tetsOnTargetSurface_array(N);
    const auto PreCalcHelper = [&]() {
        for (int i=0; i<face_on_shell.size(); i++) {
            if (t_is_removed[i])
                continue;
            for (int j=0; j<4; j++) {
                if (face_on_shell[i][j] == surfaceIdx) {
                    // this tet is on at least one face of "shellName"
                    int faceIdx = faceIdx_on_shell[i][j];
                    tetsOnTargetSurface_array[faceIdx].insert(i);

                    // later only check this candidate set for the vertices on edge0x
                    vertsOnTargetEdge_candidate_array[faceIdx].insert(TO[i][(j+1)%4]);
                    vertsOnTargetEdge_candidate_array[faceIdx].insert(TO[i][(j+2)%4]);
                    vertsOnTargetEdge_candidate_array[faceIdx].insert(TO[i][(j+3)%4]);
                }
            }
        }
    };
    PreCalcHelper();
    logger().info("Pre-calculation done: tetsOnTargetSurface & vertsOnTargetEdge");

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
        GetTetFromPrism(VO, TO, tetsOnTargetSurface_array[i], vertsOnTargetEdge_candidate_array[i], face_on_shell, surfaceIdx, prism, dualShell.V, map_VI2VO, t_is_removed,  // const input
                        T_temp, is_surface_facet_temp, face_on_shell_temp, labels_temp_vec, singularCnt);  // output
    }

    // report singularity
    logger().info("Singularity count for shell #{}: type1 = {}, type2 = {}, type3 = {}", shellName, singularCnt[0], singularCnt[1], singularCnt[2]);

    // we also want to enlarge "t_is_removed" accordingly
    // and of course no tet should be discarded
    std::vector<bool> falseVector(T_temp.size(), false);
    t_is_removed.insert( t_is_removed.end(), falseVector.begin(), falseVector.end() );

    // Update labels
    if (!args.shell_type_debug) {
        // All the tets generated here must be in the same region
        if (surfaceIdx == SURFACE_INNER)
            labels_temp = Eigen::VectorXi::Ones(T_temp.size(), 1).array() * SHELL_INNER_BOTTOM;
        else if (surfaceIdx == SURFACE_OUTER)
            labels_temp = Eigen::VectorXi::Ones(T_temp.size(), 1).array() * SHELL_TOP_OUTER;
    } else {
        labels_temp.resize(labels_temp_vec.size(), 1);
        for (int i=0; i<labels_temp_vec.size(); i++)
            labels_temp(i) = labels_temp_vec[i];
    }
}


void UpdateVertexAttributes(std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO) {

    for (int i=0; i<VO.size(); i++) {
        VO[i].conn_tets.clear();  // reset
    }

    for (int i=0; i<TO.size(); i++) {
        for (int j=0; j<4; j++) {
            VO[TO[i][j]].conn_tets.insert(i);
        }
    }

    // no need to update on_edge / on_face
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
    for (int i=0; i<F1.rows(); i++) {
        // Make sure F[i, 0] is the smallest among the three indices
        int smallIdx = 0;
        if (F1(i, 1) < F1(i, smallIdx)) smallIdx = 1;
        if (F1(i, 2) < F1(i, smallIdx)) smallIdx = 2;
        if (smallIdx != 0) {
            int t = F1(i, 0);
            F1(i, 0) = F1(i, smallIdx);
            F1(i, smallIdx) = t;
        }
    }
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
    const tetwild::Args &args,
    const DualShell_t &dualShell, 
    std::vector<tetwild::TetVertex> &VO,
    std::vector<std::array<int, 4>> &TO, 
    Eigen::VectorXi &labels,
    std::vector<std::array<int, 4>> &is_surface_facet, 
    std::vector<std::array<int, 4>> &face_on_shell, 
    std::vector<std::array<int, 4>> &faceIdx_on_shell) {

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
    Eigen::VectorXi temp;

    //// FOR shell_inner_bottom
    GenTetMeshFromShell(args, VO, TO, dualShell, face_on_shell, faceIdx_on_shell, SHELL_INNER_BOTTOM, T_temp, labels_temp, is_surface_facet_temp, face_on_shell_temp, t_is_removed);
    // concatenate new tet with old
    TO.insert( TO.end(), T_temp.begin(), T_temp.end() );
    is_surface_facet.insert( is_surface_facet.end(), is_surface_facet_temp.begin(), is_surface_facet_temp.end() );
    face_on_shell.insert( face_on_shell.end(), face_on_shell_temp.begin(), face_on_shell_temp.end() );
      temp = labels;
      labels.resize(temp.rows() + labels_temp.rows());
      labels << temp, labels_temp;
    logger().debug("GenTetMeshFromShell: SHELL_INNER_BOTTOM done");

    /// FOR shell_top_outer
    GenTetMeshFromShell(args, VO, TO, dualShell, face_on_shell, faceIdx_on_shell, SHELL_TOP_OUTER, T_temp, labels_temp, is_surface_facet_temp, face_on_shell_temp, t_is_removed);
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
    // Update "TetVertex" attributes
    UpdateVertexAttributes(VO, TO);

    logger().info("Replace with prismatic tetrahedra done");
}


void GetMeshWithPseudoTets(const DualShell_t &dualShell, const std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO, std::vector<std::array<int, 4>> &TO_with_pseudo_tets) {

    const shell_t &hallow_shell = dualShell.shell_bottom_top;
    std::map<int, int> map_VI2VO;

    TO_with_pseudo_tets = TO;

    // get the map from VI-index to VO-index
    MapIndex(VO, dualShell, dualShell.shell_bottom_top, map_VI2VO);

    for (int i=0; i<hallow_shell.size(); i++) {

        const prism_t &prism = hallow_shell[i];
        auto tetSplit = (prism[1] > prism[2]) ? TETRA_SPLIT_A : TETRA_SPLIT_B;
        if (!(prism[0] < prism[1] && prism[0] < prism[2])) {
            tetwild::log_and_throw("GetMeshWithPseudoTets: Invalid prism. prism[0] not the smallest index.");
        }

        // insert the three tets
        for (int j=0; j<3; j++) {

            std::array<int, 4> newTet;
            newTet[0] = map_VI2VO.at(prism[tetSplit[j][0]]);
            newTet[1] = map_VI2VO.at(prism[tetSplit[j][1]]);
            newTet[2] = map_VI2VO.at(prism[tetSplit[j][2]]);
            newTet[3] = map_VI2VO.at(prism[tetSplit[j][3]]);

            // degenerated due to singularity
            if (IsDegeneratedTet(VO[newTet[0]].pos, VO[newTet[1]].pos, VO[newTet[2]].pos, VO[newTet[3]].pos)) {
                // logger().warn("{} {} {} {}", newTet[0], newTet[1], newTet[2], newTet[3]);
                continue;
            }

            TO_with_pseudo_tets.push_back(newTet);
        }
    }
}


void FreezeVertices(const std::vector<std::array<int, 4>> &face_on_shell, const std::vector<std::array<int, 4>> &TO, std::vector<tetwild::TetVertex> &VO, std::vector<std::array<int, 4>> &is_surface_facet) {

    // update
    for (int i=0; i<face_on_shell.size(); i++)
        for (int j=0; j<4; j++) {

            int vIdx1 = TO[i][(j+1)%4];
            int vIdx2 = TO[i][(j+2)%4];
            int vIdx3 = TO[i][(j+3)%4];
            // Only lock vertices on SURFACE_BOTTOM and SURFACE_TOP
            if (face_on_shell[i][j] == SURFACE_BOTTOM || face_on_shell[i][j] == SURFACE_TOP) {

                // VO[vIdx].is_locked = true;
                // VO[vIdx].is_on_surface = true;

                // freeze vertices
                VO[vIdx1].is_frozen = true;
                VO[vIdx2].is_frozen = true;
                VO[vIdx3].is_frozen = true;

                // update is_surface_facet
                if (is_surface_facet[i][j] != 1 && is_surface_facet[i][j] != -1) {
                    logger().warn("is_surface_facet[{}][{}] = {}", i, j, is_surface_facet[i][j]);
                }
                VO[vIdx1].is_on_surface = true;
                VO[vIdx2].is_on_surface = true;
                VO[vIdx3].is_on_surface = true;

                // update surface_type
                if (face_on_shell[i][j] == SURFACE_BOTTOM) {
                    VO[vIdx1].surface_type = 1;
                    VO[vIdx2].surface_type = 1;
                    VO[vIdx3].surface_type = 1;
                } else {
                    VO[vIdx1].surface_type = 2;
                    VO[vIdx2].surface_type = 2;
                    VO[vIdx3].surface_type = 2;
                }
            } else if (face_on_shell[i][j] == SURFACE_INNER || face_on_shell[i][j] == SURFACE_OUTER) {
                // For other two surfaces
                VO[vIdx1].is_on_surface = false;
                VO[vIdx2].is_on_surface = false;
                VO[vIdx3].is_on_surface = false;
                is_surface_facet[i][j] = 1024;  // state.NOT_SURFACE; FIXME
            } else {
                // pass
            }
        }
}

}  // namespace tetshell
