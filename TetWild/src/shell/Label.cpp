#include <shell/Label.h>
#include <Eigen/Dense>


namespace tetshell {

void LabelTet(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI, const Eigen::MatrixXd &VO, const Eigen::MatrixXi &TO, Eigen::VectorXi &labels) {

    const int N = TO.rows();
    labels.resize(N, 1);
    
    for (int i=0; i<N; i++) {
        if (VO[FO[i, 0], 0] < -0.5)
            labels[i] = 0;
        else if (VO[FO[i, 0], 0] < 0)
            labels[i] = 1;
        else if (VO[FO[i, 0], 0] < 0.5)
            labels[i] = 2;
        else
            labels[i] = 3;
    }
}

}  // namespace tetshell
