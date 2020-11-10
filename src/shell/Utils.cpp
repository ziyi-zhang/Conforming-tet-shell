#include <shell/Shell.h>
#include <shell/Utils.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>
#include <igl/boundary_facets.h>
#include <igl/edges.h>

#include <pymesh/MshSaver.h>
#include <Eigen/Dense>
#include <unordered_map>
#include <unordered_set>


namespace tetshell {

using tetwild::Point_3;
using tetwild::logger;

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////

bool IsTetPositive(const Point_3 &p1, const Point_3 &p2, const Point_3 &p3, const Point_3 &p4) {

    CGAL::Orientation ori;
    ori = CGAL::orientation(p1, p2, p3, p4);

    return (ori == CGAL::POSITIVE);
}


bool IsTetPositive(const std::array<Point_3, 4> verts) {

    return IsTetPositive(verts[0], verts[1], verts[2], verts[3]);
}


bool MakeTetPositive(const std::vector<tetwild::TetVertex> &VO, std::array<int, 4> &T) {

    if (!IsTetPositive(VO[T[0]].pos, VO[T[1]].pos, VO[T[2]].pos, VO[T[3]].pos)) {
        int t = T[0]; T[0] = T[2]; T[2] = t;
        return true;
    }
    return false;
}


bool MakeTetPositive(const std::vector<tetwild::TetVertex> &VO, std::array<int, 4> &T, std::array<int, 4> &is_fs, std::array<int, 4> &face) {

    if (!IsTetPositive(VO[T[0]].pos, VO[T[1]].pos, VO[T[2]].pos, VO[T[3]].pos)) {
        int t;
        t = T[0]; T[0] = T[2]; T[2] = t;
        t = is_fs[0]; is_fs[0] = is_fs[2]; is_fs[2] = t;
        t = face[0]; face[0] = face[2]; face[2] = t;
        return true;
    }
    return false;
}


void UnorderedsetIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::unordered_set<int>& s) {

    if (s2.size() < s1.size()) { 
        UnorderedsetIntersection(s2, s1, s); 
        return; 
    }

    s.clear();
    s.reserve(std::min(s1.size(), s2.size()));
    for (int x : s1) {
        if (s2.count(x)) {
            s.insert(x);
        }
    }
}


void UnorderedsetIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, const std::unordered_set<int>& s3, std::unordered_set<int>& s) {

    std::unordered_set<int> s_temp;
    UnorderedsetIntersection(s1, s2, s_temp);
    UnorderedsetIntersection(s_temp, s3, s);
}


int EulerNumber(const std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO, int &V, int &F, int &E) {

    std::unordered_set<int> Vindex;  // not sure whether all vertices in VO are used
    std::unordered_set<std::pair<int, int> > Eindex;
    
    int t[4];
    for (int i=0; i<TO.size(); i++) {


        t[0] = TO[i][0];
        t[1] = TO[i][1];
        t[2] = TO[i][2];
        t[3] = TO[i][3];
        // sort 
        for (int ii=0; ii<3; ii++)
            for (int jj=ii; jj<4; jj++) {
                if (t[ii]>t[jj]) {
                    int temp = t[ii];
                    t[ii] = t[jj];
                    t[jj] = temp;
                }
            }
        // count vertices
        Vindex.insert(t[0]);
        Vindex.insert(t[1]);
        Vindex.insert(t[2]);
        Vindex.insert(t[3]);
        // count edges
        Eindex.insert(std::make_pair(t[0], t[1]));
        Eindex.insert(std::make_pair(t[0], t[2]));
        Eindex.insert(std::make_pair(t[0], t[3]));
        Eindex.insert(std::make_pair(t[1], t[2]));
        Eindex.insert(std::make_pair(t[1], t[3]));
        Eindex.insert(std::make_pair(t[2], t[3]));
    }
    V = Vindex.size();
    E = Eindex.size();

    // count faces
    Eigen::MatrixXi TO_eigen(TO.size(), 4);
    Eigen::MatrixXi F_boundary;
    for (int i=0; i<TO.size(); i++) {
        TO_eigen(i, 0) = TO[i][0];
        TO_eigen(i, 1) = TO[i][1];
        TO_eigen(i, 2) = TO[i][2];
        TO_eigen(i, 3) = TO[i][3];
    }
    igl::boundary_facets(TO_eigen, F_boundary);
    F = F_boundary.rows();

    logger().info("Euler V={} F={} E={} V+F-E={}", V, F, E, V+F-E);

    return V+F-E;
}


int EulerNumber(const std::vector<tetwild::TetVertex> &tet_vertices, const std::vector<std::array<int, 4>> &tet_indices) {

    int V, F, E;
    return EulerNumber(tet_vertices, tet_indices, V, F, E);
}


void ExtractMesh(
    const std::vector<tetwild::TetVertex> &VI, 
    const std::vector<std::array<int, 4>> &TI, 
    const Eigen::VectorXi &LI,
    const std::vector<bool> &t_is_removed,
    Eigen::MatrixXd &V_out, 
    Eigen::MatrixXi &T_out, 
    Eigen::VectorXi &L_out) {

    const int tetNum = std::count(t_is_removed.begin(), t_is_removed.end(), false);

    // v_ids is the vector of index of vertices
    std::vector<int> v_ids;
    for (int i=0; i<TI.size(); i++) {
        if (t_is_removed[i])
            continue;
        for (int j=0; j<4; j++)
            v_ids.push_back(TI[i][j]);
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

    // Prepare V, T, L
    V_out.resize(v_ids.size(), 3);
    T_out.resize(tetNum, 4);
    L_out.resize(tetNum);
    // Fill V
    for (int i = 0; i < v_ids.size(); i++) {
        for (int j = 0; j < 3; j++) {
            V_out(i, j) = CGAL::to_double(VI[v_ids[i]].pos[j]);
        }
    }
    // Fill T & L
    int cnt = 0;
    for (int i = 0; i < TI.size(); i++) {
        if (t_is_removed[i]) {
            continue;
        }
        for (int j = 0; j < 4; j++) {
            T_out(cnt, j) = map_ids[TI[i][j]];  // the index is new as defined in map_ids
        }
        L_out(cnt) = LI(i);
        cnt++;
    }

    logger().debug("Final output #VO = {}", V_out.rows());
    logger().debug("Final output #TO = {}", T_out.rows());
}


void SaveToTetMsh(const std::string fileName, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXi &L) {

    PyMesh::MshSaver mSaver(fileName, true);
    PyMesh::VectorF V_flat(V.size());
    PyMesh::VectorI T_flat(T.size());
    Eigen::MatrixXd VV = V.transpose();
    Eigen::MatrixXi TT = T.transpose();
    std::copy_n(VV.data(), V.size(), V_flat.data());
    std::copy_n(TT.data(), T.size(), T_flat.data());
    mSaver.save_mesh(V_flat, T_flat, 3, mSaver.TET);
    // mSaver.save_elem_scalar_field("min_dihedral_angle", A);
    mSaver.save_elem_scalar_field("label", L.cast<double>());

    logger().info("Result msh saved to {}", fileName);
}

}  // namespace tetshell
