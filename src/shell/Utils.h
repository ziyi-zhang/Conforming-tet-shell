#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <tetwild/CGALTypes.h>
#include <tetwild/TetmeshElements.h>
#include <pymesh/MshSaver.h>

#include <Eigen/Dense>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

void ExtractMesh(const std::vector<tetwild::TetVertex> &tet_vertices, const std::vector<std::array<int, 4>> &tet_indices, const Eigen::VectorXi &labels, const std::vector<bool> &t_is_removed, Eigen::MatrixXd &V_out, Eigen::MatrixXi &T_out, Eigen::VectorXi &L_out);
/// TODO

void SaveToTetMsh(const std::string fileName, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXi &L);

}  // namespace tetshell
