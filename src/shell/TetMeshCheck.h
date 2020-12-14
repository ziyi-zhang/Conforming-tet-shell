#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <tetwild/TetmeshElements.h>
#include <tetwild/CGALTypes.h>

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <set>
#include <limits>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

typedef struct TetMeshCheckArgs_t {
    bool positiveTet = true;  // whether all tets are positive
    bool conform = true;  // whether preserve user input
    bool vertexAttri = true;  // whether vectex attributes are correct
    bool boundary = true;  // check the boundary of region 1 and 3
    bool euler = true;  // compute the global Euler
} TetMeshCheckArgs_t;


////////////////////////////////////////////////
// Sanity check of the output tetrahedral mesh

class TetMeshCheck {
public:
    const TetMeshCheckArgs_t &args;
    const Eigen::MatrixXd &VI;
    const Eigen::MatrixXi &FI;
          std::vector<tetwild::TetVertex> &VO;  // not const b/c rational->double conversion
    const std::vector<std::array<int, 4>> &TO;
    const Eigen::VectorXi &labels;
    const std::vector<std::array<int, 4>> face_on_shell;

    bool SanityCheck(int eulerNumber = std::numeric_limits<int>::max());
    bool ConformityCheck();

    bool PositiveTetCheck();
    bool ConformCheck();
    bool VertexAttriCheck();
    bool BoundaryCheck();
    bool EulerCheck();

    void ConvertDouble();

    // Maintenance methods
    TetMeshCheck(const Eigen::MatrixXd &VI_, const Eigen::MatrixXi &FI_, std::vector<tetwild::TetVertex> &VO_, 
                 const std::vector<std::array<int, 4>> &TO_, const Eigen::VectorXi &labels_, 
                 const std::vector<std::array<int, 4>> &face_on_shell_, const TetMeshCheckArgs_t &args_) : args(args_), VI(VI_), FI(FI_), VO(VO_), TO(TO_), labels(labels_), face_on_shell(face_on_shell_) {}

    ~TetMeshCheck() {}
};

}  // namespace tetshell
