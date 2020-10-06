// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once


#include <Eigen/Dense>

#include <nlohmann/json.hpp>

namespace tetshell {

    // Json
    using json = nlohmann::json;

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

    typedef Eigen::Matrix<double, 3, 3> Matrix3;

    typedef Eigen::Matrix<double, 3, 1> Vector3;
    typedef Eigen::Matrix<double, 2, 1> Vector2;


    typedef Eigen::Matrix<int, 4, 1> Vector4i;
    typedef Eigen::Matrix<int, 3, 1> Vector3i;
    typedef Eigen::Matrix<int, 2, 1> Vector2i;
}
