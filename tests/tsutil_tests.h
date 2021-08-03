//
// Created by achains on 03.08.2021.
//

#ifndef TWIST_N_SYNC_CPP_MODULE_TSUTIL_TESTS_H
#define TWIST_N_SYNC_CPP_MODULE_TSUTIL_TESTS_H

#include "Eigen/Core"
#include "gtest/gtest.h"

#include <vector>

namespace tsutil_tests{
    // There T1, T2 are expected to be std::vector or Eigen::Vector
    template<typename T1, typename T2>
    testing::AssertionResult compareVectors(T1 const & vec1, T2 const & vec2, double eps = 0.0);

    testing::AssertionResult compareMatrices(std::vector<std::vector<double>> const & vec1, Eigen::MatrixX3d const & vec2);
}

#endif //TWIST_N_SYNC_CPP_MODULE_TSUTIL_TESTS_H
