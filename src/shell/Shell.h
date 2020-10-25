#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

typedef struct DualShell_t {

    RowMatX3d V_inner, V_bottom, V_top, V_outer;
    RowMatX3i F;
} DualShell_t;

////////////////////////////////////////////////
// Dual shell related operation

void RetrieveFourSurfaces(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI, DualShell_t &dualShell);
/// TODO


}  // namespace tetshell
