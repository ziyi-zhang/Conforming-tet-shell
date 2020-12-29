#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <tetwild/TetmeshElements.h>
#include <tetwild/CGALTypes.h>

#include <Eigen/Dense>
#include <vector>
#include <array>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

////////////////////////////////////////////////
// Estimate the optimal energy for a tetrahedron

bool EstimateOptimalEnergy(const tetwild::Point_3 &pt1, const tetwild::Point_3 &pt2, const tetwild::Point_3 &pt3, const tetwild::Point_3 &pt4, 
                           bool pt1Frozen, bool pt2Frozen, bool pt3Frozen, bool pt4Frozen, double &energy);
//bool EstimateOptimalEnergy(const std::array<double, 3> &pt1, const std::array<double, 3> &pt2, const std::array<double, 3> &pt3, const std::array<double, 3> &pt4, 
//                           bool pt1Frozen, bool pt2Frozen, bool pt3Frozen, bool pt4Frozen, double &energy);
/// TODO

void OptimizeTetWithRing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const int tetIdx, const int vertIdx);
/// TODO

}  // namespace tetshell
