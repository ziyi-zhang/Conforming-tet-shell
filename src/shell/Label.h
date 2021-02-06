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
/// return whether "point" is inside a prism
/// "tetra_split_AB" controls how to split the prism "verts"

bool PointInShell(const tetwild::Point_3 &center, const shell_t &shell, const std::vector<tetwild::Point_3> &VI);
/// return whether "center" is inside the specified "shell"
/// Brute force: VERY SLOW

void LabelTet(const tetwild::Args &args, const std::vector<tetwild::Point_3> &VI, const Eigen::MatrixXi &FI, const std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO, const std::vector<std::array<int, 4>> &face_on_shell, DualShell_t &dualShell, Eigen::VectorXi &labels);
/// Label the region of every tetrahedron in {VO, TO}
/// 
///         Region-0
///  ---- OUTER-SURFACE ---- (4)
///         Region-3
///  ----  TOP-SURFACE  ---- (3)
///         Region-2
///  ---  BOTTOM-SURFACE --- (2)
///         Region-1
///  ---  INNER-SURFACE  --- (1)
///         Region-0

void LabelInOut(const std::vector<tetwild::TetVertex> &tet_vertices, const std::vector<std::array<int, 4>> &tet_indices, const std::vector<bool> &t_is_removed, Eigen::VectorXi &labels);
/// Label tets bounded by BOTTOM_SURFACE as 1, tets bounded by TOP_SURFACE and bounding box as 2
/// Region-2 should be void when labelling
/// There may be cavities
/// [NOTE: we need this because we are not maintaining labels in the optimization stage]

}  // namespace tetshell
