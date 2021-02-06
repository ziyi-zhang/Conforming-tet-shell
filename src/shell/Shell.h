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

///         Region-0
///  ---- OUTER-SURFACE ---- (4)
///         Region-3
///  ----  TOP-SURFACE  ---- (3)
///         Region-2
///  ---  BOTTOM-SURFACE --- (2)
///         Region-1
///  ---  INNER-SURFACE  --- (1)
///         Region-0
// Terrible coding practice
// Region
#define SHELL_INNER_BOTTOM 1
#define SHELL_BOTTOM_TOP 2
#define SHELL_TOP_OUTER 3
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
/// Re-organize "VI" and "FI" to "dualShell" struct

void ReplaceWithPrismTet(const tetwild::Args &args, const DualShell_t &dualShell, std::vector<tetwild::TetVertex> &VO, std::vector<std::array<int, 4>> &TO, Eigen::VectorXi &labels, std::vector<std::array<int, 4>> &is_surface_facet, std::vector<std::array<int, 4>> &face_on_shell, std::vector<std::array<int, 4>> &faceIdx_on_shell);
/// Read in the tetrahedralized "VO", "TO" with labels
/// (1) delete all tets whose label is non-zero: all tets between inner and outer surface
/// (2) fill region-1 and region-3 with new tets split from prisms: the input surfaces are a prismatic shell, and each prism can be split into 3 tets
/// (3) update all tet labels & vertex labels. Leave region-2 empty

void GetMeshWithPseudoTets(const DualShell_t &dualShell, const std::vector<tetwild::TetVertex> &VO, const std::vector<std::array<int, 4>> &TO, std::vector<std::array<int, 4>> &TO_with_pseudo_tets);
/// Fill the empty region-2. This is not an inplace operation.
/// It is only used to calculate the Euler number for sanity check

void FreezeVertices(const std::vector<std::array<int, 4>> &face_on_shell, const std::vector<std::array<int, 4>> &TO, std::vector<tetwild::TetVertex> &VO, std::vector<std::array<int, 4>> &is_surface_facet);
/// Mark the vertices and faces on bottom and top surface as frozen
/// They cannot move when doing optimization

}  // namespace tetshell
