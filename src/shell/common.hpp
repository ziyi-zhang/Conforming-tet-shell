#ifndef PRISM_COMMON_HPP
#define PRISM_COMMON_HPP

#include <tetwild/CGALTypes.h>
#include <Eigen/Core>
#include <vector>
#include <set>
#include <map>
#include <array>
#include <tuple>
/////////////////////////////////////////////
// TypeDefs
/////////////////////////////////////////////

using RowMatd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
using RowMati = Eigen::Matrix<int, -1, -1, Eigen::RowMajor>;
using RowMat3d = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
using RowMat34i = Eigen::Matrix<size_t, 3, 4, Eigen::RowMajor>;
using RowMat34d = Eigen::Matrix<double, 3, 4, Eigen::RowMajor>;
using RowMat43d = Eigen::Matrix<double, 4, 3, Eigen::RowMajor>;
using RowMatX3d = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using RowMatX3i = Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>;
using RowMatX4i = Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajor>;
using Vec3d = Eigen::RowVector3d;
using Vec3i = std::array<int, 3>;
using Vec4i = std::array<int, 4>;

/////////////////////////////////////////////
// Prism Split A/B
/////////////////////////////////////////////

// TYPEA
// [0,3,4,5] # top
// [1,4,2,0], # bottom
// [2,5,0,4] # middle
// Boundary Facets (Outward) [0,1,4], [0,4,3], [0,3,5], [0,5,2], [0,2,1],
// [3,4,5]
//                AND [1,2,4], [2,5,4]
// Internal Facets (Upwards): [0,1,2], [0,4,2], [0,4,5], [3,4,5]
//
// TYPE B
// [0,3,4,5] # top
// [1,4,5,0] # middle
// [2,5,0,1] # bottom
// Boundary Facets (Outward) [0,1,4], [0,4,3], [0,3,5], [0,5,2], [0,2,1],
// [3,4,5] (Same As Above) Boundary Facets AND [1,2,5], [1,5,4] Internal:
// [0,1,2], [0,1,5], [0,4,5], [3,4,5]

// A reordering of the three tets, so that the first two vertices in
// each is the pillar to rely on.
constexpr std::array<std::array<size_t, 4>, 3> TETRA_SPLIT_A{
    {{0, 3, 4, 5}, {1, 4, 2, 0}, {2, 5, 0, 4}}};
// {{{0, 1, 2, 4}, {0, 4, 2, 5}, {0, 4, 5, 3}}};
constexpr std::array<std::array<size_t, 4>, 3> TETRA_SPLIT_B{
    {{0, 3, 4, 5}, {2, 5, 0, 1}, {1, 4, 5, 0}}};
// {{{0, 1, 2, 5}, {0, 1, 5, 4}, {0, 4, 5, 3}}};

constexpr std::array<std::array<size_t, 3>, 8> PRISM_BOUNDARY_A{
    {{0, 1, 4},
     {0, 4, 3},
     {0, 3, 5},
     {0, 5, 2},
     {0, 2, 1},
     {3, 4, 5}, /*the 6 afore is common for both A and B*/
     {1, 2, 4},
     {2, 5, 4}}};

constexpr std::array<std::array<size_t, 3>, 8> PRISM_BOUNDARY_B{
    {{0, 1, 4},
     {0, 4, 3},
     {0, 3, 5},
     {0, 5, 2},
     {0, 2, 1},
     {3, 4, 5}, /*the 6 afore is common for both A and B*/
     {1, 2, 5},
     {1, 5, 4}}};

constexpr std::array<std::array<size_t, 3>, 4> PRISM_RISING_FACETS_A{
    {{0, 1, 2}, {0, 4, 2}, {0, 4, 5}, {3, 4, 5}}};
constexpr std::array<std::array<size_t, 3>, 4> PRISM_RISING_FACETS_B{
    {{0, 1, 2}, {0, 1, 5}, {0, 4, 5}, {3, 4, 5}}};

inline std::tuple<bool, Vec3i, int> tetra_split_AorB(const Vec3i& abc) {
  // return a cyclic permuatation such that the minimum is the first,
  // and bool indicates whether it is well sorted. type B if sorted, A (default)
  // if not.
  bool typeA = false;
  int a = abc[0];
  int b = abc[1];
  int c = abc[2];

  // argmin
  int argmin = 0;
  if (a > b)
    argmin = b > c ? 2 : 1;
  else /*a<b*/
    argmin = a < c ? 0 : 2;

  // sort
  if (a > b) {
    std::swap(a, b);
    typeA = !typeA;
  }

  if (b > c) {
    std::swap(c, b);
    typeA = !typeA;
    if (a > b) {
      std::swap(a, b);
      typeA = !typeA;
    }
  }

  if (typeA) std::swap(b, c);

  if (typeA) assert(b > c);

  return std::tuple<bool, Vec3i, int>(typeA, Vec3i{a, b, c}, argmin);
}

const auto STANDARD_PRISM = std::array<Vec3d, 6>{
    Vec3d(0, 0, 0), Vec3d(1, 0, 0), Vec3d(0.5, sqrt(3) / 2, 0),
    Vec3d(0, 0, 1), Vec3d(1, 0, 1), Vec3d(0.5, sqrt(3) / 2, 1)};

// const auto STANDARD_PRISM = std::array<Vec3d, 6>{
//   Vec3d(0,0,0),
//   Vec3d(1,0,0),
//   Vec3d(0,1,0),
//   Vec3d(0,0,1),
//   Vec3d(1,0,1),
//   Vec3d(0,1,1)
// };

const auto SUBTRACT_TOP_ROW =
    (RowMat34d() << -1, 1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 1).finished();

// All the 12 possible tetrahedron inside a prism. When we measure quality and
// positivity, use all of them.
// The following Python scripts generates the combinations
// ```python
// import itertools
// import numpy as np
// V = np.array([[0,0,0],
//               [1,0,0],
//               [0,1,0],
//               [0,0,1],
//               [1,0,1],
//               [0,1,1]
//               ])
// for tet in set(itertools.combinations(range(6),4)):
//   tet = np.array(tet)
//   vol = np.linalg.det(V[tet[1:]]-V[tet[0]])
//   if vol == 0.0:
//     continue
//   if vol < 0:
//     tet[2:] = tet[-1], tet[-2]
//   print(tet)
// ```
// 11.27.2019, Reorder them to make mask easier for freezer
constexpr std::array<std::array<int, 4>, 12> TWELVE_TETRAS{{
    {0, 2, 3, 4},
    {0, 3, 4, 5},
    {0, 1, 5, 3},
    {0, 1, 2, 3},  // first four rely on pillar 0-3
    {1, 2, 3, 4},
    {0, 1, 2, 4},
    {0, 1, 5, 4},
    {1, 3, 4, 5},  // rely on pillar 1-4
    {1, 2, 3, 5},
    {0, 2, 5, 4},
    {0, 1, 2, 5},
    {2, 3, 4, 5}  // rely on 2-5
}};

#endif
