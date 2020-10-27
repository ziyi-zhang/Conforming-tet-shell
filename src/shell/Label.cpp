#include <shell/common.hpp>
#include <shell/Label.h>
#include <shell/Shell.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>
#include <tetwild/TetmeshElements.h>

#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <igl/copyleft/cgal/orient3D.h>
#include <Eigen/Dense>
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

/*
void GetTetFromPrism(const prism_t &prism, Eigen::MatrixXd &V, Eigen::MatrixXi &T) {

    V.resize(6, 3);
    T.resize(3, 4);
    for (int i=0; i<6; i++) {
        V.row(i) = prism[i];
    }
    T.row(0) << 0, 3, 4, 5;
    T.row(1) << 1, 4, 2, 0;
    T.row(2) << 2, 5, 0, 4;
}


void GenTetMeshFromShell(const shell_t &shell, int l, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::VectorXd &A, Eigen::VectorXi &label) {

    const int N = shell.size();
    V.resize(0, 3);
    T.resize(0, 4);

    Eigen::MatrixXd V_temp, V_old;
    Eigen::MatrixXi T_temp, T_old;
    for (int i=0; i<N; i++) {
        
        prism_t pr = shell[i];
        // shound we split the prism further?

        // get tet from prism, store in temp V & temp T
        GetTetFromPrism(shell[i], V_temp, T_temp);
        // concatenate to V T
        V_old = V;
        T_old = T;
        V.resize(V_old.rows() + V_temp.rows(), V_old.cols());
        T.resize(T_old.rows() + T_temp.rows(), T_old.cols());
        V << V_old, V_temp;
        T << T_old, T_temp + Eigen::MatrixXi::Constant(T_temp.rows(), T_temp.cols(), V_old.rows());
    }


    int num = T.rows();
    A.setZero(num, 1);
    label = Eigen::VectorXi::Constant(num, 1, l);
}


int FindIdx(const Eigen::MatrixXd &V, const Eigen::VectorXd &row) {

    for (int i=0; i<V.rows(); i++) {
        if (V(i, 0) == row(0) && V(i, 1) == row(1) && V(i, 2) == row(2))
            return i;
    }

    cerr << "prism vertex not found" << row << endl;
    return -1;
}
*/

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
        if (PointInShell(center, dualShell.shell_inner_bottom))
            labels(i) = 1;
        else if (PointInShell(center, dualShell.shell_bottom_top))
            labels(i) = 2;
        else if (PointInShell(center, dualShell.shell_top_outer))
            labels(i) = 3;
    }

    logger().info("Tet label done");
}


void ReplaceWithPrismTet(
    const DualShell_t &dualShell, 
    std::vector<tetwild::TetVertex> &VO, 
    std::vector<std::array<int, 4>> &TO, 
    Eigen::VectorXi &labels,
    std::vector<bool> &t_is_removed) {

    const int numTet = TO.size();

    // Remove labelled tets
    int numOldTet = 0;
    for (int i=0; i<numTet; i++) {
        if (labels(i) != 0) 
            t_is_removed[i] = true;
        else
            numOldTet++;
    }
    logger().debug("#tets after removing labeled shell tets = {}", numOldTet);

    // Generate new tets
/*
    // Put the new tets in VO
    Eigen::MatrixXd V_temp;
    Eigen::MatrixXi T_temp;
    Eigen::VectorXd A_temp, A_old;
    Eigen::VectorXi label_temp, label_old;
    // 1st shell
    cerr << "doing 1st shell" << endl;
    GetTetFromShell(shell_inner_bottom, 1, V_temp, T_temp, A_temp, label_temp);
    cerr << "done tet from shell" << endl;
    UnionTetMesh(VO_, TO_, V_temp, T_temp);
    cerr << "done union" << endl;
                // igl::copyleft::cgal::mesh_boolean();
    A_old = AO_;
    AO_.resize(A_old.rows() + A_temp.rows(), A_old.cols());
    AO_ << A_old, A_temp;
    label_old = label_;
    label_.resize(label_old.rows() + label_temp.rows(), label_old.cols());
    label_ << label_old, label_temp;
    // 2nd shell
    // this is user input, do not insert

    // 3rd shell
    cerr << "doing 3rd shell" << endl;
    GetTetFromShell(shell_top_outer, 3, V_temp, T_temp, A_temp, label_temp);
    UnionTetMesh(VO_, TO_, V_temp, T_temp);
    A_old = AO_;
    AO_.resize(A_old.rows() + A_temp.rows(), A_old.cols());
    AO_ << A_old, A_temp;
    label_old = label_;
    label_.resize(label_old.rows() + label_temp.rows(), label_old.cols());
    label_ << label_old, label_temp;
*/
}

}  // namespace tetshell
