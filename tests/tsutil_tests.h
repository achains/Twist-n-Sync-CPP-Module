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
    testing::AssertionResult compareVectors(T1 const & vec1, T2 const & vec2, double eps = 0.0){
        if (vec1.size() != vec2.size())
            return testing::AssertionFailure() << "vec1.size() != vec2.size() <<>> "
            << vec1.size() << " != " << vec2.size();
        for (Eigen::Index i = 0; i < vec1.size(); ++i)
            if (std::abs(vec1[i] - vec2[i]) > eps)
                return testing::AssertionFailure() << "Elements at index " << i << " are not equal with eps=" << eps;
        return testing::AssertionSuccess();
    }
}

#endif //TWIST_N_SYNC_CPP_MODULE_TSUTIL_TESTS_H
