#include <shell/common.hpp>
#include <shell/Label.h>

#include <igl/copyleft/cgal/orient3D.h>
#include <Eigen/Dense>
#include <array>


namespace tetshell {

using namespace std;

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
    /*
    cerr << "doing 2nd shell" << endl;
    GetTetFromShell(shell_bottom_top, 2, V_temp, T_temp, A_temp, label_temp);
    UnionTetMesh(VO_, TO_, V_temp, T_temp);
    A_old = AO_;
    AO_.resize(A_old.rows() + A_temp.rows(), A_old.cols());
    AO_ << A_old, A_temp;
    label_old = label_;
    label_.resize(label_old.rows() + label_temp.rows(), label_old.cols());
    label_ << label_old, label_temp;
    */
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

}  // namespace tetshell
