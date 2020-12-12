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


    logger().warn("Input shell has problematic singularites.");
    return true;
}

}  // namespace tetshell
