#include <shell/Shell.h>
#include <shell/Utils.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>

#include <pymesh/MshSaver.h>
#include <Eigen/Dense>
#include <unordered_map>


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


void MakeTetPositive(const std::vector<tetwild::TetVertex> &VO, std::array<int, 4> &T) {

    if (!IsTetPositive(VO[T[0]].pos, VO[T[1]].pos, VO[T[2]].pos, VO[T[3]].pos)) {
        int t = T[0]; T[0] = T[2]; T[2] = t;
    }
}


void MakeTetPositive(const std::vector<tetwild::TetVertex> &VO, std::array<int, 4> &T, std::array<int, 4> &is_fs, std::array<int, 4> &face) {

    if (!IsTetPositive(VO[T[0]].pos, VO[T[1]].pos, VO[T[2]].pos, VO[T[3]].pos)) {
        int t;
        t = T[0]; T[0] = T[2]; T[2] = t;
        t = is_fs[0]; is_fs[0] = is_fs[2]; is_fs[2] = t;
        t = face[0]; face[0] = face[2]; face[2] = t;
    }
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

    logger().debug("final output #v = {}", V_out.rows());
    logger().debug("final output #t = {}", T_out.rows());
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
    std::cout << "Result msh saved to " << fileName << std::endl;
}

}  // namespace tetshell