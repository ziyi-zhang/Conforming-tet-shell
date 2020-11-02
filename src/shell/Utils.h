#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <tetwild/CGALTypes.h>
#include <tetwild/TetmeshElements.h>
#include <pymesh/MshSaver.h>

#include <Eigen/Dense>
#include <vector>
#include <unordered_set>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

bool IsTetPositive(const tetwild::Point_3 &p1, const tetwild::Point_3 &p2, const tetwild::Point_3 &p3, const tetwild::Point_3 &p4);
bool IsTetPositive(const std::array<tetwild::Point_3, 4> verts);
/// TODO

bool MakeTetPositive(const std::vector<tetwild::TetVertex> &VO, std::array<int, 4> &T);
bool MakeTetPositive(const std::vector<tetwild::TetVertex> &VO, std::array<int, 4> &T, std::array<int, 4> &is_fs, std::array<int, 4> &face);
/// TODO

void UnorderedsetIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::unordered_set<int>& s);
/// TODO

void ExtractMesh(const std::vector<tetwild::TetVertex> &tet_vertices, const std::vector<std::array<int, 4>> &tet_indices, const Eigen::VectorXi &labels, const std::vector<bool> &t_is_removed, Eigen::MatrixXd &V_out, Eigen::MatrixXi &T_out, Eigen::VectorXi &L_out);
/// TODO

void SaveToTetMsh(const std::string fileName, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXi &L);
/// TODO

}  // namespace tetshell
