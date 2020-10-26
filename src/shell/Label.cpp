#include <shell/common.hpp>
#include <shell/Label.h>
#include <shell/Shell.h>
#include <tetwild/CGALTypes.h>

#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <igl/copyleft/cgal/orient3D.h>
#include <Eigen/Dense>
#include <array>


namespace tetshell {

using namespace std;
using tetwild::Point_3;

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


void GetTetFromShell(const shell_t &shell, int l, Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::VectorXd &A, Eigen::VectorXi &label) {

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


void UnionTetMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::MatrixXd &V_temp, Eigen::MatrixXi &T_temp) {

    const int N = V_temp.rows();
    Eigen::VectorXi lookup(N);

    for (int i=0; i<N; i++) {
        // where is V_temp(i) in V?
        lookup(i) = FindIdx(V, V_temp.row(i));
        if (lookup(i) == -1) {
            V.conservativeResize(V.rows()+1, V.cols());
            V.row(V.rows()-1) = V_temp.row(i);  // insert a new vertex
            lookup(i) = V.rows()-1;
        }
    }

    for (int i=0; i<T_temp.rows(); i++) {
        
        T_temp(i, 0) = lookup(T_temp(i, 0));
        T_temp(i, 1) = lookup(T_temp(i, 1));
        T_temp(i, 2) = lookup(T_temp(i, 2));
        T_temp(i, 3) = lookup(T_temp(i, 3));
    }

    Eigen::MatrixXi T_old = T;
    T.resize(T_old.rows() + T_temp.rows(), T_old.cols());
    T << T_old, T_temp;
}
*/

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LabelTet(
    const std::vector<Point_3> &VI, 
    const Eigen::MatrixXi &FI, 
    const std::vector<tetwild::TetVertex> &VO, 
    const std::vector<std::array<int, 4>> &TO,
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
}

/*
void ReplaceWithPrismTet(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI, Eigen::MatrixXd &VO, Eigen::MatrixXi &TO, Eigen::VectorXd &AO, Eigen::VectorXi &labels) {

    const int M = VI.rows() / 4;
    const int Ntri = FI.rows() / 4;
    const int Ntet = TO.rows();

    cerr << "starting to remove" << endl;
    // Remove all labelled tets
    std::vector<bool> mask(Ntet, true);
    int numOldTet = 0, count = 0;
    for (int i=0; i<Ntet; i++) {
        if (labels(i) != 0) 
            mask[i] = false;
        else
            numOldTet++;
    }
    cerr << "numOldTet " << numOldTet << endl;
    Eigen::MatrixXd VO_ = VO;
    Eigen::MatrixXi TO_(numOldTet, 4);
    Eigen::VectorXd AO_(numOldTet, 1);
    Eigen::VectorXi label_ = Eigen::VectorXi::Constant(numOldTet, 1, 0);
    for (int i=0; i<Ntet; i++) {
        if (mask[i]) {
            // survived
            TO_.row(count) = TO.row(i);
            AO_(count) = AO(i);
            count++;
        }
    }
    cerr << "remove done" << endl;

    ///// TRASH code...
    // generate new tets
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


    // replace inplace
    VO = VO_;
    TO = TO_;
    AO = AO_;
    labels = label_;
}
*/

}  // namespace tetshell
