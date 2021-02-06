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
// Estimate the optimal (smallest) energy of tetrahedron {pt1, pt2, pt3, pt4} when {pt1Frozen+pt2Frozen+pt3Frozen+pt4Frozen} vertices are frozen
// If all four vertices are locked, the energy is a constant
// NOTE: we are ignoring the other tetrahedra when estimating the energy, so this is a quite rough estimation.

void OptimizeTetWithRing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const int tetIdx, const int vertIdx);
/// DEBUG PURPOSE

}  // namespace tetshell
