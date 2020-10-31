#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <tetwild/CGALTypes.h>

#include <Eigen/Dense>
#include <vector>
#include <set>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

// FIXME: move this to State as a local struct
// Region
#define SHELL_INNER_BOTTOM 1
#define SHELL_BOTTOM_TOP 2
#define SHELL_TOP_OUTER 3
// FIXME: move this to State as a local struct
// Surface
#define NOT_SUR 0
#define SURFACE_INNER 1
#define SURFACE_BOTTOM 2
#define SURFACE_TOP 3
#define SURFACE_OUTER 4


typedef std::array<int, 6> prism_t;  // stores the vectex index A1,A2,A3+B1,B2,B3
typedef std::vector<prism_t> shell_t;

typedef struct DualShell_t {

    std::vector<tetwild::Point_3> V;  // V_inner, V_bottom, V_top, V_outer concatenated
    RowMatX3i F;  // for each of the four V, not all their concatenation

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

bool PointInShell(const tetwild::Point_3 &center, const shell_t &shell, const std::vector<tetwild::Point_3> &VI);
/// TODO

// void SaveMsh();


}  // namespace tetshell
