#include <shell/Utils.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>
#include <tetwild/TetmeshElements.h>

#include <Eigen/Dense>
#include <algorithm>


namespace tetshell {

using tetwild::Point_3;
using tetwild::Point_3f;
using tetwild::logger;

namespace {

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


void tetwild_comformalAMIPSJacobian_new(const double * T, double *result_0) {
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
    double helper_1 = helper_0[1];
    double helper_2 = helper_0[10];
    double helper_3 = helper_1 - helper_2;
    double helper_4 = helper_0[0];
    double helper_5 = helper_0[3];
    double helper_6 = helper_0[9];
    double helper_7 = 0.577350269189626*helper_4 - 1.15470053837925*helper_5 + 0.577350269189626*helper_6;
    double helper_8 = helper_0[2];
    double helper_9 = 0.408248290463863*helper_8;
    double helper_10 = helper_0[5];
    double helper_11 = 0.408248290463863*helper_10;
    double helper_12 = helper_0[8];
    double helper_13 = 1.22474487139159*helper_12;
    double helper_14 = helper_0[11];
    double helper_15 = 0.408248290463863*helper_14;
    double helper_16 = helper_11 - helper_13 + helper_15 + helper_9;
    double helper_17 = 0.577350269189626*helper_8;
    double helper_18 = 1.15470053837925*helper_10;
    double helper_19 = 0.577350269189626*helper_14;
    double helper_20 = helper_17 - helper_18 + helper_19;
    double helper_21 = helper_0[6];
    double helper_22 = -1.22474487139159*helper_21 + 0.408248290463863*helper_4 + 0.408248290463863*helper_5 + 0.408248290463863*helper_6;
    double helper_23 = helper_16*helper_7 - helper_20*helper_22;
    double helper_24 = -helper_14 + helper_8;
    double helper_25 = 0.408248290463863*helper_1;
    double helper_26 = helper_0[4];
    double helper_27 = 0.408248290463863*helper_26;
    double helper_28 = helper_0[7];
    double helper_29 = 1.22474487139159*helper_28;
    double helper_30 = 0.408248290463863*helper_2;
    double helper_31 = helper_25 + helper_27 - helper_29 + helper_30;
    double helper_32 = helper_31*helper_7;
    double helper_33 = 0.577350269189626*helper_1;
    double helper_34 = 1.15470053837925*helper_26;
    double helper_35 = 0.577350269189626*helper_2;
    double helper_36 = helper_33 - helper_34 + helper_35;
    double helper_37 = helper_22*helper_36;
    double helper_38 = helper_4 - helper_6;
    double helper_39 = helper_23*helper_3 - helper_24*(helper_32 - helper_37) - helper_38*(helper_16*helper_36 - helper_20*helper_31);
    double helper_40 = pow(pow(helper_39, 2), -0.333333333333333);
    double helper_41 = 0.707106781186548*helper_10 - 0.707106781186548*helper_12;
    double helper_42 = 0.707106781186548*helper_26 - 0.707106781186548*helper_28;
    double helper_43 = 0.5*helper_21 + 0.5*helper_5;
    double helper_44 = 0.5*helper_26 + 0.5*helper_28;
    double helper_45 = 0.5*helper_10 + 0.5*helper_12;
    double helper_46 = 0.666666666666667*(helper_1*(-1.5*helper_1 + 0.5*helper_2 + helper_44) + helper_10*(-1.5*helper_10 + 0.5*helper_12 + 0.5*helper_14 + 0.5*helper_8) + helper_12*(0.5*helper_10 - 1.5*helper_12 + 0.5*helper_14 + 0.5*helper_8) + helper_14*(-1.5*helper_14 + helper_45 + 0.5*helper_8) + helper_2*(0.5*helper_1 - 1.5*helper_2 + helper_44) + helper_21*(-1.5*helper_21 + 0.5*helper_4 + 0.5*helper_5 + 0.5*helper_6) + helper_26*(0.5*helper_1 + 0.5*helper_2 - 1.5*helper_26 + 0.5*helper_28) + helper_28*(0.5*helper_1 + 0.5*helper_2 + 0.5*helper_26 - 1.5*helper_28) + helper_4*(-1.5*helper_4 + helper_43 + 0.5*helper_6) + helper_5*(0.5*helper_21 + 0.5*helper_4 - 1.5*helper_5 + 0.5*helper_6) + helper_6*(0.5*helper_4 + helper_43 - 1.5*helper_6) + helper_8*(0.5*helper_14 + helper_45 - 1.5*helper_8))/helper_39;
    double helper_47 = -0.707106781186548*helper_21 + 0.707106781186548*helper_5;
    result_0[0] = -helper_40*(1.0*helper_21 - 3.0*helper_4 + helper_46*(helper_41*(-helper_1 + helper_2) - helper_42*(helper_14 - helper_8) - (-helper_17 + helper_18 - helper_19)*(-helper_25 - helper_27 + helper_29 - helper_30) + (-helper_33 + helper_34 - helper_35)*(-helper_11 + helper_13 - helper_15 - helper_9)) + 1.0*helper_5 + 1.0*helper_6);
    result_0[1] = helper_40*(3.0*helper_1 - 1.0*helper_2 - 1.0*helper_26 - 1.0*helper_28 + helper_46*(helper_23 + helper_24*helper_47 - helper_38*helper_41));
    result_0[2] = helper_40*(-1.0*helper_10 - 1.0*helper_12 - 1.0*helper_14 + helper_46*(-helper_3*helper_47 - helper_32 + helper_37 + helper_38*helper_42) + 3.0*helper_8);
}


void tetwild_comformalAMIPSHessian_new(const double * T, double *result_0) {
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
    double helper_3 = helper_1 - helper_2;
    double helper_4 = helper_0[0];
    double helper_5 = 0.577350269189626*helper_4;
    double helper_6 = helper_0[3];
    double helper_7 = 1.15470053837925*helper_6;
    double helper_8 = helper_0[9];
    double helper_9 = 0.577350269189626*helper_8;
    double helper_10 = helper_5 - helper_7 + helper_9;
    double helper_11 = helper_0[1];
    double helper_12 = 0.408248290463863*helper_11;
    double helper_13 = helper_0[4];
    double helper_14 = 0.408248290463863*helper_13;
    double helper_15 = helper_0[7];
    double helper_16 = 1.22474487139159*helper_15;
    double helper_17 = helper_0[10];
    double helper_18 = 0.408248290463863*helper_17;
    double helper_19 = helper_12 + helper_14 - helper_16 + helper_18;
    double helper_20 = helper_10*helper_19;
    double helper_21 = 0.577350269189626*helper_11;
    double helper_22 = 1.15470053837925*helper_13;
    double helper_23 = 0.577350269189626*helper_17;
    double helper_24 = helper_21 - helper_22 + helper_23;
    double helper_25 = 0.408248290463863*helper_4;
    double helper_26 = 0.408248290463863*helper_6;
    double helper_27 = helper_0[6];
    double helper_28 = 1.22474487139159*helper_27;
    double helper_29 = 0.408248290463863*helper_8;
    double helper_30 = helper_25 + helper_26 - helper_28 + helper_29;
    double helper_31 = helper_24*helper_30;
    double helper_32 = helper_3*(helper_20 - helper_31);
    double helper_33 = helper_4 - helper_8;
    double helper_34 = 0.408248290463863*helper_1;
    double helper_35 = helper_0[5];
    double helper_36 = 0.408248290463863*helper_35;
    double helper_37 = helper_0[8];
    double helper_38 = 1.22474487139159*helper_37;
    double helper_39 = 0.408248290463863*helper_2;
    double helper_40 = helper_34 + helper_36 - helper_38 + helper_39;
    double helper_41 = helper_24*helper_40;
    double helper_42 = 0.577350269189626*helper_1;
    double helper_43 = 1.15470053837925*helper_35;
    double helper_44 = 0.577350269189626*helper_2;
    double helper_45 = helper_42 - helper_43 + helper_44;
    double helper_46 = helper_19*helper_45;
    double helper_47 = helper_41 - helper_46;
    double helper_48 = helper_33*helper_47;
    double helper_49 = helper_11 - helper_17;
    double helper_50 = helper_10*helper_40;
    double helper_51 = helper_30*helper_45;
    double helper_52 = helper_50 - helper_51;
    double helper_53 = helper_49*helper_52;
    double helper_54 = helper_32 + helper_48 - helper_53;
    double helper_55 = pow(helper_54, 2);
    double helper_56 = pow(helper_55, -0.333333333333333);
    double helper_57 = 1.0*helper_27 - 3.0*helper_4 + 1.0*helper_6 + 1.0*helper_8;
    double helper_58 = 0.707106781186548*helper_13;
    double helper_59 = 0.707106781186548*helper_15;
    double helper_60 = helper_58 - helper_59;
    double helper_61 = helper_3*helper_60;
    double helper_62 = 0.707106781186548*helper_35 - 0.707106781186548*helper_37;
    double helper_63 = helper_49*helper_62;
    double helper_64 = helper_47 + helper_61 - helper_63;
    double helper_65 = 1.33333333333333/helper_54;
    double helper_66 = 1.0/helper_55;
    double helper_67 = 0.5*helper_27 + 0.5*helper_6;
    double helper_68 = -1.5*helper_4 + helper_67 + 0.5*helper_8;
    double helper_69 = 0.5*helper_4 + helper_67 - 1.5*helper_8;
    double helper_70 = -1.5*helper_27 + 0.5*helper_4 + 0.5*helper_6 + 0.5*helper_8;
    double helper_71 = 0.5*helper_27 + 0.5*helper_4 - 1.5*helper_6 + 0.5*helper_8;
    double helper_72 = 0.5*helper_13 + 0.5*helper_15;
    double helper_73 = -1.5*helper_11 + 0.5*helper_17 + helper_72;
    double helper_74 = 0.5*helper_11 - 1.5*helper_17 + helper_72;
    double helper_75 = 0.5*helper_11 + 0.5*helper_13 - 1.5*helper_15 + 0.5*helper_17;
    double helper_76 = 0.5*helper_11 - 1.5*helper_13 + 0.5*helper_15 + 0.5*helper_17;
    double helper_77 = 0.5*helper_35 + 0.5*helper_37;
    double helper_78 = -1.5*helper_1 + 0.5*helper_2 + helper_77;
    double helper_79 = 0.5*helper_1 - 1.5*helper_2 + helper_77;
    double helper_80 = 0.5*helper_1 + 0.5*helper_2 + 0.5*helper_35 - 1.5*helper_37;
    double helper_81 = 0.5*helper_1 + 0.5*helper_2 - 1.5*helper_35 + 0.5*helper_37;
    double helper_82 = helper_1*helper_78 + helper_11*helper_73 + helper_13*helper_76 + helper_15*helper_75 + helper_17*helper_74 + helper_2*helper_79 + helper_27*helper_70 + helper_35*helper_81 + helper_37*helper_80 + helper_4*helper_68 + helper_6*helper_71 + helper_69*helper_8;
    double helper_83 = 0.444444444444444*helper_66*helper_82;
    double helper_84 = helper_66*helper_82;
    double helper_85 = -helper_32 - helper_48 + helper_53;
    double helper_86 = 1.0/helper_85;
    double helper_87 = helper_86*pow(pow(helper_85, 2), -0.333333333333333);
    double helper_88 = 0.707106781186548*helper_6;
    double helper_89 = 0.707106781186548*helper_27;
    double helper_90 = helper_88 - helper_89;
    double helper_91 = 0.666666666666667*helper_10*helper_40 + 0.666666666666667*helper_3*helper_90 - 0.666666666666667*helper_30*helper_45 - 0.666666666666667*helper_33*helper_62;
    double helper_92 = -3.0*helper_11 + 1.0*helper_13 + 1.0*helper_15 + 1.0*helper_17;
    double helper_93 = -helper_11 + helper_17;
    double helper_94 = -helper_1 + helper_2;
    double helper_95 = -helper_21 + helper_22 - helper_23;
    double helper_96 = -helper_34 - helper_36 + helper_38 - helper_39;
    double helper_97 = -helper_42 + helper_43 - helper_44;
    double helper_98 = -helper_12 - helper_14 + helper_16 - helper_18;
    double helper_99 = -0.666666666666667*helper_60*helper_94 + 0.666666666666667*helper_62*helper_93 + 0.666666666666667*helper_95*helper_96 - 0.666666666666667*helper_97*helper_98;
    double helper_100 = helper_3*helper_90;
    double helper_101 = helper_33*helper_62;
    double helper_102 = helper_100 - helper_101 + helper_52;
    double helper_103 = -helper_60*helper_94 + helper_62*helper_93 + helper_95*helper_96 - helper_97*helper_98;
    double helper_104 = 0.444444444444444*helper_102*helper_103*helper_82*helper_86 + helper_57*helper_91 - helper_92*helper_99;
    double helper_105 = 1.85037170770859e-17*helper_1*helper_78 + 1.85037170770859e-17*helper_11*helper_73 + 1.85037170770859e-17*helper_13*helper_76 + 1.85037170770859e-17*helper_15*helper_75 + 1.85037170770859e-17*helper_17*helper_74 + 1.85037170770859e-17*helper_2*helper_79 + 1.85037170770859e-17*helper_27*helper_70 + 1.85037170770859e-17*helper_35*helper_81 + 1.85037170770859e-17*helper_37*helper_80 + 1.85037170770859e-17*helper_4*helper_68 + 1.85037170770859e-17*helper_6*helper_71 + 1.85037170770859e-17*helper_69*helper_8;
    double helper_106 = helper_64*helper_82*helper_86;
    double helper_107 = -0.666666666666667*helper_10*helper_19 + 0.666666666666667*helper_24*helper_30 + 0.666666666666667*helper_33*helper_60 - 0.666666666666667*helper_49*helper_90;
    double helper_108 = -3.0*helper_1 + 1.0*helper_2 + 1.0*helper_35 + 1.0*helper_37;
    double helper_109 = -helper_20 + helper_31 + helper_33*helper_60 - helper_49*helper_90;
    double helper_110 = 0.444444444444444*helper_109*helper_82*helper_86;
    double helper_111 = helper_103*helper_110 + helper_107*helper_57 - helper_108*helper_99;
    double helper_112 = -helper_4 + helper_8;
    double helper_113 = -helper_88 + helper_89;
    double helper_114 = -helper_5 + helper_7 - helper_9;
    double helper_115 = -helper_25 - helper_26 + helper_28 - helper_29;
    double helper_116 = helper_82*helper_86*(helper_112*helper_62 + helper_113*helper_94 + helper_114*helper_96 - helper_115*helper_97);
    double helper_117 = -helper_100 + helper_101 - helper_50 + helper_51;
    double helper_118 = -helper_102*helper_110 + helper_107*helper_92 + helper_108*helper_91;
    double helper_119 = helper_82*helper_86*(helper_112*(-helper_58 + helper_59) - helper_113*helper_93 - helper_114*helper_98 + helper_115*helper_95);
    result_0[0] = helper_56*(helper_57*helper_64*helper_65 - pow(helper_64, 2)*helper_83 + 0.666666666666667*helper_64*helper_84*(-helper_41 + helper_46 - helper_61 + helper_63) + 3.0);
    result_0[1] = helper_87*(helper_104 - helper_105*helper_35 + helper_106*helper_91);
    result_0[2] = helper_87*(helper_106*helper_107 + helper_111);
    result_0[3] = helper_87*(helper_104 + helper_116*helper_99);
    result_0[4] = helper_56*(-pow(helper_117, 2)*helper_83 + helper_117*helper_65*helper_92 + helper_117*helper_84*helper_91 + 3.0);
    result_0[5] = helper_87*(-helper_105*helper_6 - helper_107*helper_116 + helper_118);
    result_0[6] = helper_87*(-helper_105*helper_13 + helper_111 + helper_119*helper_99);
    result_0[7] = helper_87*(helper_118 - helper_119*helper_91);
    result_0[8] = helper_56*(-helper_108*helper_109*helper_65 - 1.11111111111111*pow(helper_109, 2)*helper_84 + 3.0);
}


bool NewtonsUpdate(
    const std::array<Point_3f, 4> &pts,
    const int idx,
    double& energy, 
    Eigen::Vector3d &J, 
    Eigen::Matrix3d &H, 
    Eigen::Vector3d &X_old) {
// Adapted from VertexSmoother

    // initialize and update x_old
    for (int i=0; i<3; i++) {
        J(i) = 0.0;
        for (int j=0; j<3; j++)
            H(i, j) = 0.0;
        X_old(i) = pts[idx][i];
    }

    // stores the coordinates of 4 3D-vertices to pos
    std::array<double, 12> pos;
    for (int j=0; j<4; j++) {
        for (int k=0; k<3; k++) {
            pos[j*3+k] = pts[(idx+j) % 4][k];
        }
    }

    double J_1[3];
    double H_1[9];
    tetwild_comformalAMIPSJacobian_new(pos.data(), J_1);
    tetwild_comformalAMIPSHessian_new(pos.data(), H_1);
    energy = tetwild_comformalAMIPSEnergy_new(pos.data());

    for (int j=0; j<3; j++) {
        J(j) += J_1[j];
        H(j, 0) += H_1[j * 3 + 0];
        H(j, 1) += H_1[j * 3 + 1];
        H(j, 2) += H_1[j * 3 + 2];
    }

    if (std::isinf(energy)) {
        energy = 1e50;
    }
    if (std::isnan(energy)) {
        return false;
    }
    if (energy <= 0) {
        return false;
    }
    if (!J.allFinite()) {
        return false;
    }
    if (!H.allFinite()) {
        return false;
    }

    return true;
}


bool NewtonsUpdateBatch(
    const std::vector<Point_3f> &pts, 
    const std::vector<std::array<int, 4> > &tets, 
    const int targetTetIdx,  // only used to track the interested tet's energy
    const int targetVertIdx, 
    double &energy, 
    Eigen::Vector3d &J, 
    Eigen::Matrix3d &H, 
    Eigen::Vector3d &X_old) {

    // initialize
    for (int i=0; i<3; i++) {
        J(i) = 0.0;
        for (int j=0; j<3; j++)
            H(i, j) = 0.0;
        X_old(i) = pts[targetVertIdx][i];
    }

    for (int i=0; i<tets.size(); i++) {

        // stores the coordinates of 4 3D-vertices to pos
        std::array<double, 12> pos;
        int start = -1;  // this is necessary to get the correct sign of J
        for (int j=0; j<4; j++) {
            if (tets[i][j] == targetVertIdx) {
                start = j;
                break;
            }
        }
        for (int j=0; j<4; j++) {
            for (int k=0; k<3; k++) {
                pos[j*3+k] = pts[tets[i][(start+j) % 4]][k];
            }
        }

        double J_1[3];
        double H_1[9];
        tetwild_comformalAMIPSJacobian_new(pos.data(), J_1);
        tetwild_comformalAMIPSHessian_new(pos.data(), H_1);
        if (i == targetTetIdx)
            energy = tetwild_comformalAMIPSEnergy_new(pos.data());

        for (int j=0; j<3; j++) {
            J(j) += J_1[j];
            H(j, 0) += H_1[j * 3 + 0];
            H(j, 1) += H_1[j * 3 + 1];
            H(j, 2) += H_1[j * 3 + 2];
        }
    }

    if (std::isinf(energy)) {
        logger().debug("NewtonsUpdateBatch: std::isinf(energy)");
        energy = 1e50;
    }
    if (std::isnan(energy)) {
        logger().debug("NewtonsUpdateBatch: std::isnan(energy)");
        return false;
    }
    if (energy <= 0) {
        logger().debug("NewtonsUpdateBatch: energy <= 0");
        return false;
    }
    if (!J.allFinite()) {
        logger().debug("NewtonsUpdateBatch: !J.allFinite()");
        return false;
    }
    if (!H.allFinite()) {
        logger().debug("NewtonsUpdateBatch: !H.allFinite()");
        return false;
    }

    return true;
}


bool OneStepNewtonOptim(std::array<Point_3, 4> &pts_rational, int idx, double &energy) {

    const int stepMax = 100;  // line search iteration
    Eigen::Vector3d J;
    Eigen::Matrix3d H;
    Eigen::Vector3d X_old;
    double oldEnergy, newEnergy;

    // to double
    std::array<Point_3f, 4> pts;
    for (int i=0; i<4; i++)
        pts[i] = Point_3f(CGAL::to_double(pts_rational[i][0]), CGAL::to_double(pts_rational[i][1]), CGAL::to_double(pts_rational[i][2]));

    // jacobian and hessian
    if (NewtonsUpdate(pts, idx, oldEnergy, J, H, X_old) == false)
        return false;

    // Initialize line search
    Point_3f old_pf = pts[idx];
    Point_3 old_p = pts_rational[idx];
    bool step_taken = false;

    // solve linear system
    Eigen::Vector3d delta_X = H.colPivHouseholderQr().solve(J);

    // whether the solution is valid
    bool newtons = true;
    if (!delta_X.allFinite()) {
        newtons = false;
    } else if (!J.isApprox(H * delta_X, 1e-10)) {
        newtons = false;
    }

    // use newton's or gradient descent
    Eigen::Vector3d optimDirection;
    if (newtons) {
        optimDirection = delta_X.stableNormalized();
    } else {
        optimDirection = J.stableNormalized();
    }

    // naive line search
    double a = 1.0;
    for (int step=0; step<stepMax; step++) {

        // try to update location to be [ X_old - a * optimDirection ]  
        pts[idx] = Point_3f(X_old(0) - a*optimDirection(0), X_old(1) - a*optimDirection(1), X_old(2) - a*optimDirection(2));
        pts_rational[idx] = Point_3(X_old(0) - a*optimDirection(0), X_old(1) - a*optimDirection(1), X_old(2) - a*optimDirection(2));

        // check flipping
        if (!IsTetPositive(pts_rational[0], pts_rational[1], pts_rational[2], pts_rational[3])) {
            // abort change
            pts[idx] = old_pf;
            pts_rational[idx] = old_p;
            a /= 2.0;
            continue;
        }

        // check energy
        // stores the coordinates of 4 3D-vertices to pos
        std::array<double, 12> pos;
        for (int j=0; j<4; j++) {
            for (int k=0; k<3; k++) {
                pos[j*3+k] = pts[j][k];
            }
        }
        newEnergy = tetwild_comformalAMIPSEnergy_new(pos.data());
        if (newEnergy >= oldEnergy || std::isinf(newEnergy) || std::isnan(newEnergy)) {
            // abort change
            pts[idx] = old_pf;
            pts_rational[idx] = old_p;
            a /= 2.0;
            continue;
        }

        step_taken = true;
        energy = newEnergy;
        break;
    }  // for (int it=0; it<MAX_IT; it++)

    if (J.norm() < 1e-9 && std::abs(newEnergy - oldEnergy) < 1e-5)
        step_taken = false;  // break if there is hardly anything left can be done

    return step_taken;
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////

bool EstimateOptimalEnergy(const tetwild::Point_3 &pt1, const tetwild::Point_3 &pt2, const tetwild::Point_3 &pt3, const tetwild::Point_3 &pt4, 
                           bool pt1Frozen, bool pt2Frozen, bool pt3Frozen, bool pt4Frozen, double &energy) {

    // Valid input?
    if (!IsTetPositive(pt1, pt2, pt3, pt4)) {
        logger().warn("EstimateOptimalEnergy: Input tetrahedron not positive.");
        return false;
    }

    // Newton's
    double old_energy = energy;  // debug
    std::array<Point_3, 4> pts_rational{pt1, pt2, pt3, pt4};
    std::array<bool, 4> frozenArray{pt1Frozen, pt2Frozen, pt3Frozen, pt4Frozen};
    const int roundMax = 300;
    bool moved = false;
    int roundCnt;
    for (roundCnt=0; roundCnt<roundMax; roundCnt++) {

        bool optimizedInRound = false;
        for (int vIdx=0; vIdx<4; vIdx++) {

            if (frozenArray[vIdx]) continue;  // the only constrain is the locked vertices
            if (OneStepNewtonOptim(pts_rational, vIdx, energy)) {
                optimizedInRound = true;
                moved = true;
            }
        }

        if (!optimizedInRound) {
            break;
        }
    }

    int cntFrozen = int(pt1Frozen == true) + int(pt2Frozen == true) + int(pt3Frozen == true) + int(pt4Frozen == true);
    // logger().critical("{} round={} frozen={} energy={} old_e={}", moved, roundCnt, cntFrozen, energy, old_energy);

    return moved;
}


void OptimizeTetWithRing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const int targetTetIdx, const int vertIdx) {

    // delete unused V & create rational copy
    std::vector<Point_3f> pts;
    std::vector<Point_3>  pts_rational;
    std::vector<std::array<int, 4> > tets;
    int targetVertIdx = -1;
    for (int i=0; i<T.rows(); i++) {

        std::array<int, 4> newTetIdx;

        // for the j-th tet
        for (int j=0; j<4; j++) {
            Point_3f pt{V(T(i, j), 0), V(T(i, j), 1), V(T(i, j), 2)};
            auto it = find(pts.begin(), pts.end(), pt);
            if (it == pts.end()) {
                // a new vertex
                pts.push_back(pt);
                pts_rational.push_back( Point_3(pt[0], pt[1], pt[2]) );
                newTetIdx[j] = pts.size() - 1;
            } else {
                // already exist
                newTetIdx[j] = it - pts.begin();
            }
            if (i==targetTetIdx && j==vertIdx) targetVertIdx = newTetIdx[j];
        }
        // insert the updated T
        tets.push_back(newTetIdx);
    }


    // Start optimization
    const int itMax = 20;  // newton's iteration
    const int stepMax = 100;  // line search iteration
    Eigen::Vector3d J;
    Eigen::Matrix3d H;
    Eigen::Vector3d X_old;
    double oldEnergy = -1.0, newEnergy;
    for (int it=0; it<itMax; it++) {

        // update jacobian, hessian and X_old
        if (!NewtonsUpdateBatch(pts, tets, targetTetIdx, targetVertIdx, oldEnergy, J, H, X_old)) {
            logger().warn("NewtonsUpdateBatch returned false with it={}. Now exit optimization.", it);
            break;
        }
        logger().debug("@Start of it={} energy={}", it, oldEnergy);

        // initialize
        Point_3f old_pf = pts[targetVertIdx];
        Point_3 old_p = pts_rational[targetVertIdx];
        bool step_taken = false;

        // solve linear system
        Eigen::Vector3d delta_X = H.colPivHouseholderQr().solve(J);
        logger().debug("it={}  delta_X=[{} {} {}]", it, delta_X[0], delta_X[1], delta_X[2]);  // DEBUG PURPOSE

        // whether the solution is valid
        bool newtons = true;
        if (!delta_X.allFinite()) {
            newtons = false;
            logger().debug("delta_X.allFinite() not satisfied. Will use gradient instead.");
        } else if (!J.isApprox(H * delta_X, 1e-10)) {
            newtons = false;
            logger().debug("J.isApprox(H * delta_X, 1e-10) not satisfied. Will use gradient instead.");
        }

        // use newton's or gradient descent
        Eigen::Vector3d optimDirection;
        if (newtons) {
            optimDirection = delta_X.stableNormalized();
        } else {
            optimDirection = J.stableNormalized();
            logger().debug("J=[{} {} {}]", J[0], J[1], J[2]);
        }

        // naive line search
        double a = 1.0;
        for (int step=0; step<stepMax; step++) {

            // try to update location to be [ X0 - a * optimDirection ]
            pts[targetVertIdx] = Point_3f(X_old(0) - a*optimDirection(0), X_old(1) - a*optimDirection(1), X_old(2) - a*optimDirection(2));
            pts_rational[targetVertIdx] = Point_3(X_old(0) - a*optimDirection(0), X_old(1) - a*optimDirection(1), X_old(2) - a*optimDirection(2));

            // check flipping
            if (!AreTetsPositive(pts_rational, tets)) {
                // abort change
                logger().debug("   it={} step={} a={} rejected due to flipped tet", it, step, a);
                pts[targetVertIdx] = old_pf;
                pts_rational[targetVertIdx] = old_p;
                a /= 2.0;
                continue;
            }

            // check energy
            // stores the coordinates of 4 3D-vertices to pos
            std::array<double, 12> pos;
            for (int j=0; j<4; j++) {
                for (int k=0; k<3; k++) {
                    pos[j*3+k] = pts[tets[targetTetIdx][j]][k];
                }
            }
            newEnergy = tetwild_comformalAMIPSEnergy_new(pos.data());
            if (newEnergy >= oldEnergy || std::isinf(newEnergy) || std::isnan(newEnergy)) {
                // abort change
                logger().debug("   it={} step={} a={} newE={} oldE={} rejected due to energy", it, step, a, newEnergy, oldEnergy);
                pts[targetVertIdx] = old_pf;
                pts_rational[targetVertIdx] = old_p;
                a /= 2.0;
                continue;
            }

            step_taken = true;
            break;
        }  // for (int step=0; step<stepMax; step++)  line search loop
        if (!step_taken)
            logger().debug("Line search failed to find a valid step length @ it={}", it);

        // if (std::abs(newEnergy - oldEnergy) < 1e-8)
        //    step_taken = false;
        if (!step_taken) {
            break;
        }
    }  // for (int it=0; it<itMax; it++)  optimization iter loop
    logger().debug("Final energy = {}", newEnergy);
}

}  // namespace tetshell
