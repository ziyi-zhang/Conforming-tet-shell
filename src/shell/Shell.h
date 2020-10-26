#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <tetwild/CGALTypes.h>

#include <Eigen/Dense>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

typedef std::array<tetwild::Point_3, 6> prism_t;
typedef std::vector<prism_t> shell_t;

typedef struct DualShell_t {

    std::vector<tetwild::Point_3> V_inner, V_bottom, V_top, V_outer;
    RowMatX3i F;

    shell_t shell_inner_bottom, shell_bottom_top, shell_top_outer;
} DualShell_t;

////////////////////////////////////////////////
// Dual shell related operation

void GenDualShell(const std::vector<tetwild::Point_3> &VI, const Eigen::MatrixXi &FI, DualShell_t &dualShell);
/// TODO

bool point_in_tetrahedron(const tetwild::Point_3& point, const tetwild::Point_3& T0, const tetwild::Point_3& T1, const tetwild::Point_3& T2, const tetwild::Point_3& T3);
/// return whether "point" is inside a tetrahedron formed by T0, T1, T2, T3
/// especially, "on_boundary" does not count as inside
/// note: this function does not require T0, T1, T2, T3 to form a positive tetrahedron

bool point_in_prism(const tetwild::Point_3& point, bool tetra_split_AB, const std::array<tetwild::Point_3, 6>& verts);
/// TODO

bool PointInShell(const tetwild::Point_3 &center, const shell_t &shell);
/// TODO

void SaveMsh();


}  // namespace tetshell
