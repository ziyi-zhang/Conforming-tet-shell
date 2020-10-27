#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <shell/Shell.h>
#include <tetwild/TetmeshElements.h>
#include <tetwild/CGALTypes.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

////////////////////////////////////////////////
// Label tetrahedra based on their position

void LabelTet(const std::vector<tetwild::Point_3> &VI, const Eigen::MatrixXi &FI, const std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO, const std::vector<std::array<int, 4>> &face_on_shell, DualShell_t &dualShell, Eigen::VectorXi &labels);
/// TODO

void ConstructShellFromTri(const Eigen::MatrixXd &V_bottom, const Eigen::MatrixXd &V_top, const Eigen::MatrixXi &F, shell_t &shell);
/// TODO

void ReplaceWithPrismTet(const DualShell_t &dualShell, std::vector<tetwild::TetVertex> &VO, std::vector<std::array<int, 4>> &TO, Eigen::VectorXi &labels, std::vector<bool> &t_is_removed);
/// TODO

}  // namespace tetshell
