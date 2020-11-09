#include <shell/TetMeshCheck.h>
#include <shell/common.hpp>
#include <shell/Shell.h>
#include <shell/Utils.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>

#include <Eigen/Dense>


namespace tetshell {

using namespace std;
using tetwild::logger;

namespace {
}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////

bool TetMeshCheck::SanityCheck() {

    logger().info("======================================");
    logger().info("Tetrahedral mesh sanity check results:");
    bool result = true;

    if (args.positiveTet) {
        bool res = PositiveTetCheck();
        result = result && res;
    }

    
    logger().info("Tetrahedral mesh sanity check done.");
    logger().info("======================================");
    return result;
}


bool TetMeshCheck::PositiveTetCheck() {

    logger().debug(">>> PositiveTetCheck >>>");
    bool result = true;

    for (auto it=TO.begin(); it!=TO.end(); it++) {
        if (!IsTetPositive(VO[it->at(0)].pos, VO[it->at(1)].pos, VO[it->at(2)].pos, VO[it->at(3)].pos)) {
            logger().warn("Flipped or degenerated tetrahedron. #TO = {}", it-TO.begin());
            result = false;
        }
    }

    logger().debug("<<< PositiveTetCheck <<<");
    return result;
}

}  // namespace tetshell
