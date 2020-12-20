#include <shell/ShellCheck.h>
#include <tetwild/Logger.h>
#include <igl/boundary_facets.h>

#include <Eigen/Dense>


namespace tetshell {

using namespace std;
using tetwild::logger;

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////

bool ShellCheck::SanityCheck() {

    bool result = true;

    if (args.boundary) {
        bool res = BoundaryCheck();
        result = result && res;
    }

    if (args.singularity) {
        bool res = SingularityCheck();
        result = result && res;
    }

    if (args.Findex) {
        bool res = FindexCheck();
        result = result && res;
    }

    logger().info("Input shell sanity check done.");
    logger().info("======================================");
    return result;
}


bool ShellCheck::BoundaryCheck() {

    Eigen::MatrixXi F;
    
    igl::boundary_facets(FI, F);
    if (F.rows() > 0) {
        logger().warn("Input shell not closed. Boundary edge #={}", F.rows());
        return false;
    }

    return true;
}


bool ShellCheck::SingularityCheck() {

    int Nv = VI.rows() / 4;
    bool singularityZone = true;
    for (int i=0; i<Nv; i++) {
        if (singularityZone) {
            if (SamePoint(i, i+Nv) && SamePoint(i, i+Nv*2) && SamePoint(i, i+Nv*3)) {
                continue;
            } else {
                singularityZone = false;
                continue;
            }
        }

        // if find any collision, abort
        if (SamePoint(i, i+Nv) || SamePoint(i, i+Nv*2) || SamePoint(i, i+Nv*3) || 
            SamePoint(i+Nv, i+Nv*2) || SamePoint(i+Nv, i+Nv*3) || SamePoint(i+Nv*2, i+Nv*3)) {
                logger().warn("Input shell has problematic singularites among [{}, {}, {}, {}]", i, i+Nv, i+Nv*2, i+Nv*3);
                return false;
            }
    }

    return true;
}


bool ShellCheck::FindexCheck() {

    // first index of a triangle must be the smallest
    for (int i=0; i<FI.rows(); i++)
        if (FI(i, 0) >= FI(i, 1) || FI(i, 0) >= FI(i, 2)) {
            logger().warn("Input connectivity matrix first index not the smallest F[{}, :] = [{}, {}, {}]", i, FI(i, 0), FI(i, 1), FI(i, 2));
            return false;
        }

    // V, F size should match
    int Nf = FI.rows() / 4;
    const Eigen::MatrixXi &FI_one_shell = FI.block(0, 0, Nf, 3);
    int maxF_index = FI_one_shell.maxCoeff();
    if (maxF_index+1 != VI.rows() / 4) {
        logger().warn("Input V, F not valid. F.max()={}, V.rows()={}", maxF_index, VI.rows());
        return false;
    }

    return true;
}


bool ShellCheck::SamePoint(int x, int y) {

    if (VI(x, 0) == VI(y, 0) &&
        VI(x, 1) == VI(y, 1) &&
        VI(x, 2) == VI(y, 2))
        return true;
    else
        return false;
}

}  // namespace tetshell
