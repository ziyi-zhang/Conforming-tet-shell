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

} TetMeshCheckArgs_t;


////////////////////////////////////////////////
// Sanity check of the output tetrahedral mesh

class TetMeshCheck {
public:
    const TetMeshCheckArgs_t &args;
    const Eigen::MatrixXd &VI;
    const Eigen::MatrixXi &FI;
    const std::vector<tetwild::TetVertex> &VO;
    const std::vector<std::array<int, 4>> &TO;

    bool SanityCheck();
    bool PositiveTetCheck();

    // Maintenance methods
    TetMeshCheck(const Eigen::MatrixXd &VI_, const Eigen::MatrixXi &FI_, const std::vector<tetwild::TetVertex> &VO_, 
                 const std::vector<std::array<int, 4>> &TO_, const TetMeshCheckArgs_t &args_) : VI(VI_), FI(FI_), VO(VO_), TO(TO_), args(args_) {}

    ~TetMeshCheck() {}
};

}  // namespace tetshell
