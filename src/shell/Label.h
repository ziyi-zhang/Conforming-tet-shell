#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <tetwild/TetmeshElements.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

////////////////////////////////////////////////
// Label tetrahedra based on their position

void LabelTet(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI, const std::vector<tetwild::TetVertex> &VO, const RowMatX4i &TO, Eigen::VectorXi &labels);
/// TODO

void ConstructShellFromTri(const Eigen::MatrixXd &V_bottom, const Eigen::MatrixXd &V_top, const Eigen::MatrixXi &F, shell_t &shell);
/// TODO

void ReplaceWithPrismTet(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI, Eigen::MatrixXd &VO, Eigen::MatrixXi &TO, Eigen::VectorXd &AO, Eigen::VectorXi &labels);
/// TODO

}  // namespace tetshell
