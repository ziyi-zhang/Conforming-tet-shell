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
bool AreTetsPositive(const std::vector<tetwild::Point_3> &V, const std::vector<std::array<int, 4> > &tets);
/// TODO

bool MakeTetPositive(const std::vector<tetwild::TetVertex> &VO, std::array<int, 4> &T);
bool MakeTetPositive(const std::vector<tetwild::TetVertex> &VO, std::array<int, 4> &T, std::array<int, 4> &is_fs, std::array<int, 4> &face);
/// TODO

void UnorderedsetIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::unordered_set<int>& s);
void UnorderedsetIntersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, const std::unordered_set<int>& s3, std::unordered_set<int>& s);
/// TODO

int EulerNumber(const std::vector<std::array<int, 4>> &tet_indices, int &V, int &F, int &E, const std::string &info);
int EulerNumber(const std::vector<std::array<int, 4>> &tet_indices, const std::string &info);
/// TODO

bool SameUnorderedTriangle(const Eigen::Matrix<double, 3, 3> &triA, const Eigen::Matrix<double, 3, 3> &triB, const double epsilon);
/// TODO

void ExtractMesh(const tetwild::Args &args, const std::vector<tetwild::TetVertex> &tet_vertices, const std::vector<std::array<int, 4>> &tet_indices, const std::vector<tetwild::TetQuality> &tetQuality, const Eigen::VectorXi &labels, const std::vector<bool> &t_is_removed, Eigen::MatrixXd &V_out, Eigen::MatrixXi &T_out, Eigen::VectorXd &A_out, Eigen::VectorXi &L_out);
/// TODO

void SaveToTetMsh(const std::string fileName, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXd &A, const Eigen::VectorXi &L);
/// TODO

void PrintPoints(const tetwild::Point_3& pt1, const tetwild::Point_3& pt2, const tetwild::Point_3& pt3, const tetwild::Point_3& pt4);
/// TODO

}  // namespace tetshell
