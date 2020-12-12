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

typedef struct ShellCheckArgs_t {
    bool boundary = true;  // whether the four shells are closed (judging from F only)
    bool singularity = true;  // 
} ShellCheckArgs_t;


////////////////////////////////////////////////
// Sanity check of the input dual shell

class ShellCheck {
public:
    const ShellCheckArgs_t &args;
    const Eigen::MatrixXd &VI;
    const Eigen::MatrixXi &FI;

    bool SanityCheck();

    bool BoundaryCheck();
    bool SingularityCheck();


    // Maintenance methods
    ShellCheck(const Eigen::MatrixXd &VI_, const Eigen::MatrixXi &FI_, const ShellCheckArgs_t &args_) : args(args_), VI(VI_), FI(FI_) {}
    ~ShellCheck() {}
};

}  // namespace tetshell
