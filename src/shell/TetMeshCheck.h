#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <shell/common.hpp>
#include <tetwild/TetmeshElements.h>
#include <tetwild/CGALTypes.h>

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <set>
////////////////////////////////////////////////////////////////////////////////

namespace tetshell {

typedef struct TetMeshCheckArgs_t {
    bool positiveTet = true;  // whether all tets are positive
    bool conform = true;  // whether preserve user input
    bool vertexAttri = true;  // whether vectex attributes are correct
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

    bool SanityCheck();
    bool PositiveTetCheck();
    bool ConformCheck();
    bool VertexAttriCheck();

    void ConvertDouble();

    // Maintenance methods
    TetMeshCheck(const Eigen::MatrixXd &VI_, const Eigen::MatrixXi &FI_, std::vector<tetwild::TetVertex> &VO_, 
                 const std::vector<std::array<int, 4>> &TO_, const TetMeshCheckArgs_t &args_) : VI(VI_), FI(FI_), VO(VO_), args(args_), TO(TO_) {}

    ~TetMeshCheck() {}
};

}  // namespace tetshell
