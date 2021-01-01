// This file compares the result with TetGen
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readPLY.h>
#include <igl/boundary_facets.h>

#include <string>
#include <set>
#include <pymesh/MshSaver.h>


double tetwild_comformalAMIPSEnergy_new(const double * T) {
    double helper_0[12];
    helper_0[0] = T[0];
    helper_0[1] = T[1];
    helper_0[2] = T[2];
    helper_0[3] = T[3];
    helper_0[4] = T[4];
    helper_0[5] = T[5];
    helper_0[6] = T[6];
    helper_0[7] = T[7];
    helper_0[8] = T[8];
    helper_0[9] = T[9];
    helper_0[10] = T[10];
    helper_0[11] = T[11];
    double helper_1 = helper_0[2];
    double helper_2 = helper_0[11];
    double helper_3 = helper_0[0];
    double helper_4 = helper_0[3];
    double helper_5 = helper_0[9];
    double helper_6 = 0.577350269189626 * helper_3 - 1.15470053837925 * helper_4 + 0.577350269189626 * helper_5;
    double helper_7 = helper_0[1];
    double helper_8 = helper_0[4];
    double helper_9 = helper_0[7];
    double helper_10 = helper_0[10];
    double helper_11 = 0.408248290463863 * helper_10 + 0.408248290463863 * helper_7 + 0.408248290463863 * helper_8 -
                       1.22474487139159 * helper_9;
    double helper_12 = 0.577350269189626 * helper_10 + 0.577350269189626 * helper_7 - 1.15470053837925 * helper_8;
    double helper_13 = helper_0[6];
    double helper_14 = -1.22474487139159 * helper_13 + 0.408248290463863 * helper_3 + 0.408248290463863 * helper_4 +
                       0.408248290463863 * helper_5;
    double helper_15 = helper_0[5];
    double helper_16 = helper_0[8];
    double helper_17 = 0.408248290463863 * helper_1 + 0.408248290463863 * helper_15 - 1.22474487139159 * helper_16 +
                       0.408248290463863 * helper_2;
    double helper_18 = 0.577350269189626 * helper_1 - 1.15470053837925 * helper_15 + 0.577350269189626 * helper_2;
    double helper_19 = 0.5 * helper_13 + 0.5 * helper_4;
    double helper_20 = 0.5 * helper_8 + 0.5 * helper_9;
    double helper_21 = 0.5 * helper_15 + 0.5 * helper_16;
    return -(helper_1 * (-1.5 * helper_1 + 0.5 * helper_2 + helper_21) +
             helper_10 * (-1.5 * helper_10 + helper_20 + 0.5 * helper_7) +
             helper_13 * (-1.5 * helper_13 + 0.5 * helper_3 + 0.5 * helper_4 + 0.5 * helper_5) +
             helper_15 * (0.5 * helper_1 - 1.5 * helper_15 + 0.5 * helper_16 + 0.5 * helper_2) +
             helper_16 * (0.5 * helper_1 + 0.5 * helper_15 - 1.5 * helper_16 + 0.5 * helper_2) +
             helper_2 * (0.5 * helper_1 - 1.5 * helper_2 + helper_21) +
             helper_3 * (helper_19 - 1.5 * helper_3 + 0.5 * helper_5) +
             helper_4 * (0.5 * helper_13 + 0.5 * helper_3 - 1.5 * helper_4 + 0.5 * helper_5) +
             helper_5 * (helper_19 + 0.5 * helper_3 - 1.5 * helper_5) +
             helper_7 * (0.5 * helper_10 + helper_20 - 1.5 * helper_7) +
             helper_8 * (0.5 * helper_10 + 0.5 * helper_7 - 1.5 * helper_8 + 0.5 * helper_9) +
             helper_9 * (0.5 * helper_10 + 0.5 * helper_7 + 0.5 * helper_8 - 1.5 * helper_9)) *
           pow(pow((helper_1 - helper_2) * (helper_11 * helper_6 - helper_12 * helper_14) -
                   (-helper_10 + helper_7) * (-helper_14 * helper_18 + helper_17 * helper_6) +
                   (helper_3 - helper_5) * (-helper_11 * helper_18 + helper_12 * helper_17), 2), -0.333333333333333);
}


double GetAMIPSEnergy(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &T,
    const int tetIdx) {

    // stores the coordinates of 4 3D-vertices to pos
    std::array<double, 12> pos;
    for (int j=0; j<4; j++) {
        for (int k=0; k<3; k++) {
            pos[j*3+k] = V(T(tetIdx, j), k);
        }
    }

    return tetwild_comformalAMIPSEnergy_new(pos.data());
}


void AddBbox(Eigen::MatrixXd &Vin, Eigen::MatrixXi &Fin) {

    // TetShell: Enlarge the bbox to allow more freedom in optimization
    double pos_gap[3] = {0, 0, 0};
    double neg_gap[3] = {0, 0, 0};
    double pmin[3] = {1e10, 1e10, 1e10};
    double pmax[3] = {-1e10, -1e10, -1e10};

    // basic bbox
    for (int i=0; i<Vin.rows(); i++) {
        for (int j=0; j<3; j++) {
            if (Vin(i, j) < pmin[j]) pmin[j] = Vin(i, j);
            if (Vin(i, j) > pmax[j]) pmax[j] = Vin(i, j);
        }
    }

    // larger bbox
    for (int i=0; i<Fin.rows(); i++)
        for (int j=0; j<3; j++) {

            int vert1 = Fin(i, (j+1) % 3);
            int vert2 = Fin(i, (j+2) % 3);
            double len = (Vin.row(vert1) - Vin.row(vert2)).norm();

            for (int k=0; k<3; k++) {
                // x, y, z
                
                if (Vin(vert1, k) + len - pmax[k] > pos_gap[k]) pos_gap[k] = Vin(vert1, k) + len - pmax[k];
                if (Vin(vert2, k) + len - pmax[k] > pos_gap[k]) pos_gap[k] = Vin(vert2, k) + len - pmax[k];

                if (pmin[k] - (Vin(vert1, k) - len) > neg_gap[k]) neg_gap[k] = pmin[k] - (Vin(vert1, k) - len);
                if (pmin[k] - (Vin(vert2, k) - len) > neg_gap[k]) neg_gap[k] = pmin[k] - (Vin(vert2, k) - len);
            }
        }

    printf("Bbox x_neg_gap=%.3f y_neg_gap=%.3f z_neg_gap=%.3f\n", neg_gap[0], neg_gap[1], neg_gap[2]);
    printf("Bbox x_pos_gap=%.3f y_pos_gap=%.3f z_pos_gap=%.3f\n", pos_gap[0], pos_gap[1], pos_gap[2]);

    // update pmin & pmax
    for (int i=0; i<3; i++) {
        pmin[i] -= neg_gap[i];
        pmax[i] += pos_gap[i];
    }

    // insert 8 vertices into V & 12 faces to F
    int oldRowsV = Vin.rows();
    Vin.conservativeResize(Vin.rows()+8, Vin.cols());
    Vin.row(oldRowsV  ) << pmin[0], pmin[1], pmin[2];
    Vin.row(oldRowsV+1) << pmax[0], pmin[1], pmin[2];
    Vin.row(oldRowsV+2) << pmax[0], pmin[1], pmax[2];
    Vin.row(oldRowsV+3) << pmin[0], pmin[1], pmax[2];
    Vin.row(oldRowsV+4) << pmin[0], pmax[1], pmin[2];
    Vin.row(oldRowsV+5) << pmax[0], pmax[1], pmin[2];
    Vin.row(oldRowsV+6) << pmax[0], pmax[1], pmax[2];
    Vin.row(oldRowsV+7) << pmin[0], pmax[1], pmax[2];
    /*
    int oldRowsF = Fin.rows();
    Fin.conservativeResize(Fin.rows()+12, Fin.cols());
    Fin.row(oldRowsF   ) << 4, 7, 8;
    Fin.row(oldRowsF+1 ) << 4, 7, 3;
    Fin.row(oldRowsF+2 ) << 5, 7, 8;
    Fin.row(oldRowsF+3 ) << 5, 7, 6;
    Fin.row(oldRowsF+4 ) << 7, 2, 3;
    Fin.row(oldRowsF+5 ) << 7, 2, 6;
    Fin.row(oldRowsF+6 ) << 1, 3, 2;
    Fin.row(oldRowsF+7 ) << 1, 3, 4;
    Fin.row(oldRowsF+8 ) << 1, 8, 4;
    Fin.row(oldRowsF+9 ) << 1, 8, 5;
    Fin.row(oldRowsF+10) << 1, 6, 2;
    Fin.row(oldRowsF+11) << 1, 6, 5;
    for (int i=oldRowsF; i<=oldRowsF+11; i++) {
        for (int j=0; j<3; j++) {
            Fin(i, j) += oldRowsV - 1;
        }
    }
    */
}


int main(int argc, char *argv[]) {

    Eigen::MatrixXd Vin, V;
    Eigen::MatrixXi Fin, F;

    if (argc > 1) {
        std::string filePath(argv[1]);
        std::cout << filePath << std::endl;

        igl::readPLY(filePath, V, F);
    }
    if (V.size() == 0 || F.size() == 0) {
        std::cout << "File read error" << std::endl;
        return 0;
    }

    // process Fin
    /*
    int facesPerSurface = F.rows();
    int Nv = Vin.rows() / 4;
    Fin.resize(facesPerSurface * 2, 3);
    Eigen::MatrixXi bottom = F.array() + Nv;
    Eigen::MatrixXi top = F.array() + Nv * 2;
    Fin << bottom, top;
    */
    int facesPerSurface = F.rows();
    int Nv = V.rows() / 4;
    Fin = F;
    // Vin = V.block(Nv, 0, Nv, 3);
    Vin = V;

    // add bounding box
    // AddBbox(Vin, Fin);

    // call tetgen
    Eigen::MatrixXd Vout;
    Eigen::MatrixXi Tout;
    Eigen::MatrixXi Fout;
    const std::string switches = "pYT1e-15MqO5V";
    int returnCode = 99;
    returnCode = igl::copyleft::tetgen::tetrahedralize(Vin, Fin, switches, Vout, Tout, Fout);
    std::cout << "switches = " << switches << std::endl;
    std::cout << "tetgen returnCode = " << returnCode << std::endl;

    // AMIPS stats
    double maxE = 0.0;
    double minE = 1000.0;
    double sumE = 0.0;
    for (int i=0; i<Tout.rows(); i++) {
        double e = GetAMIPSEnergy(Vout, Tout, i);
        if (e > maxE) maxE = e;
        if (e < minE) minE = e;
        sumE += e;
    }
    printf("maxE  = %.1f\n", maxE);
    printf("minE  = %.1f\n", minE);
    printf("meanE = %.1f\n", sumE / double(Tout.rows()));

    // conformity check
    Eigen::MatrixXi F_boundary;
    Eigen::VectorXi F_index, tet_local_idx;
    igl::boundary_facets(Tout, F_boundary, F_index, tet_local_idx);
    if (F_boundary.rows() != Fin.rows()) {
        printf("Warning: F_boundary.rows() = %ld <-> Fin.rows() = %ld\n", F_boundary.rows(), Fin.rows());
    }
    std::set<int> boundaryVerts;
    for (int i=0; i<F_boundary.rows(); i++) {
        for (int j=0; j<3; j++) {
            boundaryVerts.insert(F_boundary(i, j));
        }
    }
    if (boundaryVerts.size() != Vin.rows()) {
        printf("Warning: boundaryVerts.size() = %ld <-> Vin.rows() = %ld\n", boundaryVerts.size(), Vin.rows());
    }
    // vertex location
    const int N = boundaryVerts.size();
    std::vector<bool> matched(N, false);
    for (auto it=boundaryVerts.begin(); it!=boundaryVerts.end(); it++) {
        // for the i-th vertex in Vout
        int vout_id = *it;
        for (int k=0; k<N; k++) {
            if (matched[k]) continue;
            if (Vout(vout_id, 0)==Vin(k, 0) && Vout(vout_id, 1)==Vin(k, 1) && Vout(vout_id, 2)==Vin(k, 2)) {
                matched[k] = true;
                break;
            }
        }
    }
    int matchedNum = 0;
    for (int i=0; i<matched.size(); i++)
        matchedNum += int(matched[i]);
    if (matchedNum != Vin.rows()) {
        printf("Warning: matchedNum = %d <-> Vin.rows() = %ld\n", matchedNum, Vin.rows());
    }
    printf("Check done\n");

    // export
    PyMesh::MshSaver mSaver("./tetgen_result.msh", true);
    PyMesh::VectorF V_flat(Vout.size());
    PyMesh::VectorI T_flat(Tout.size());
    Eigen::MatrixXd VV = Vout.transpose();
    Eigen::MatrixXi TT = Tout.transpose();
    std::copy_n(VV.data(), Vout.size(), V_flat.data());
    std::copy_n(TT.data(), Tout.size(), T_flat.data());
    mSaver.save_mesh(V_flat, T_flat, 3, mSaver.TET);

    return 0;
}
