#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <shell/Shell.h>
#include <tetwild/TetmeshElements.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Args.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

////////////////////////////////////////////////
// Label tetrahedra based on their position

bool point_in_tetrahedron(const tetwild::Point_3& point, const tetwild::Point_3& T0, const tetwild::Point_3& T1, const tetwild::Point_3& T2, const tetwild::Point_3& T3);
/// return whether "point" is inside a tetrahedron formed by T0, T1, T2, T3
/// especially, "on_boundary" does not count as inside
/// note: this function does not require T0, T1, T2, T3 to form a positive tetrahedron

bool point_in_prism(const tetwild::Point_3& point, bool tetra_split_AB, const std::array<tetwild::Point_3, 6>& verts);
/// TODO

bool PointInShell(const tetwild::Point_3 &center, const shell_t &shell, const std::vector<tetwild::Point_3> &VI);
/// TODO

void LabelTet(const tetwild::Args &args, const std::vector<tetwild::Point_3> &VI, const Eigen::MatrixXi &FI, const std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO, const std::vector<std::array<int, 4>> &face_on_shell, DualShell_t &dualShell, Eigen::VectorXi &labels);
/// TODO

}  // namespace tetshell
