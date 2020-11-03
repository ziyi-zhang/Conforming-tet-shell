#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <tetwild/TetmeshElements.h>
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

void ReplaceWithPrismTet(const DualShell_t &dualShell, std::vector<tetwild::TetVertex> &VO, std::vector<std::array<int, 4>> &TO, Eigen::VectorXi &labels, std::vector<std::array<int, 4>> &is_surface_facet, std::vector<std::array<int, 4>> &face_on_shell);
/// TODO

// void SaveMsh();


}  // namespace tetshell
