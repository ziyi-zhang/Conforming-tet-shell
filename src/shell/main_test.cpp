#include <shell/OptimalEnergy.h>
#include <igl/readMSH.h>

#include <Eigen/Dense>

using namespace tetshell;

int main() {

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;

    igl::readMSH("/home/zhongshi/public/shell_4layer/62286.stl.h5_.msh", V, T);
    T.resize(28, 4);
    T << 801,   794,   800, 11128,
       921,   795,   865, 11128,
         800, 11128,   864,   799,
         801,   872, 11128,   930,
         864, 11128,   872,   893,
         872,   891, 11128,   912,
         800,   801, 11128,   799,
         864,   794,   872, 11128,
         864, 11128,   915,   928,
         794,   801,   872, 11128,
       11128,   915,   896,   864,
         896, 11128,   893,   891,
       11128,   928,   864,   865,
         891,   872, 11128,   893,
         896,   891,   915, 11128,
       11128,   865,   864,   799,
       11128,   795,   799,   801,
         928,   912,   930, 11128,
         930,   872, 11128,   912,
         928,   915,   912, 11128,
         912,   891, 11128,   915,
       11128,   794,   800,   864,
         928,   930,   865, 11128,
       11128,   896,   893,   864,
         801, 11128,   921,   930,
         801,   795,   921, 11128,
         930,   921,   865, 11128, 
       11128,   795,   865,   799;

    // OptimizeTetWithRing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const int tetIdx, const int vertIdx);
    OptimizeTetWithRing(V, T, 22, 3);
}
