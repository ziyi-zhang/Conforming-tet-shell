#include <shell/common.hpp>
#include <shell/ShellCheck.h>
#include <shell/Utils.h>
#include <tetwild/Logger.h>
#include <tetwild/CGALTypes.h>
#include <igl/boundary_facets.h>

#include <Eigen/Dense>


namespace tetshell {

using namespace std;
using tetwild::logger;
using tetwild::Point_3;

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

    if (args.positiveTet) {
        bool res = PositiveTetCheck();
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


bool PrismPositveTets(const std::vector<Point_3> &VI_cgal, const Eigen::Vector3i &tri1, const Eigen::Vector3i &tri2) {

    auto tets = tri1(1)>tri1(2) ? TETRA_SPLIT_A : TETRA_SPLIT_B;
    Eigen::VectorXi verts(6, 1);
    verts << tri1, tri2;

    for (int i = 0; i < 3; i++) {
        if (!IsTetPositive(VI_cgal[verts(tets[i][0])], VI_cgal[verts(tets[i][1])], VI_cgal[verts(tets[i][2])], VI_cgal[verts(tets[i][3])]))
            return false;
    }
    return true;
}


bool ShellCheck::PositiveTetCheck() {

    // Convert VI to CGAL rational
    std::vector<Point_3> VI_cgal;
    for (int i=0; i<VI.rows(); i++) {
        VI_cgal.push_back(Point_3(VI(i, 0), VI(i, 1), VI(i, 2)));
    }

    // check positive tet between two adjacent surfaces
    int Nf = FI.rows() / 4;
    for (int i=0; i<Nf; i++) {
        if (!PrismPositveTets(VI_cgal, FI.row(i), FI.row(i+Nf))) {
            logger().warn("Input prism not consisted of positive tets: INNER-BOTTOM with row = {}", i);
            return false;
        }
        if (!PrismPositveTets(VI_cgal, FI.row(i+Nf), FI.row(i+Nf*2))) {
            logger().warn("Input prism not consisted of positive tets: BOTTOM-TOP with row = {}", i);
            return false;
        }
        if (!PrismPositveTets(VI_cgal, FI.row(i+Nf*2), FI.row(i+Nf*3))) {
            logger().warn("Input prism not consisted of positive tets: TOP_OUTER with row = {}", i);
            return false;
        }
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
