// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once

// To set the parameters related
#include <tetshell/Types.hpp>

#include <array>
#include <vector>

#include <geogram/mesh/mesh.h>

namespace tetshell {

class Parameters {

public:
    std::string log_path    = "";
    std::string input_path  = "";
    std::string output_path = "";

    bool is_quiet  = false;
    int  log_level = 0;//2;

    unsigned int num_threads = std::numeric_limits<unsigned int>::max();
};
}  // namespace tetshell
