//
// Created by achains on 03.08.2021.
//


#include "util/TSUtil.h"
#include "tsutil_tests.h"

namespace tsutil_tests{
    template <typename T1, typename T2>
    testing::AssertionResult compareVectors(T1 const & vec1, T2 const & vec2){

        if (vec1.size() != vec2.size()) return testing::AssertionFailure() << "vec1.size() != vec2.size()";

        for (Eigen::Index i = 0; i < vec1.size(); ++i)
            if (vec1[i] != vec2[i]) return testing::AssertionFailure() << "Elements at " << i << " are not equal";
        return testing::AssertionSuccess();
    }

    TEST(VectorTransform, VectorToEigen){
        std::vector<double> vec1{1., 2., 3., 4., 5., 6.};
        Eigen::VectorXd vec2 = tsutil::vectorToEigVectorXd(vec1);

        EXPECT_TRUE(compareVectors(vec1, vec2));
    }

    // Beware the definition of tsutil::vectorToEigMatrixX3d
    // Function should only work for matrices with size (N, 3)
    TEST(VectorTransform, MatrixToEigen){
        std::vector<std::vector<double>> vec1{{3., 2., 1.}, {1., 2., 3.}};
        Eigen::MatrixX3d vec2 = tsutil::vectorToEigMatrixX3d(vec1);

        EXPECT_EQ(vec1.size(), vec2.rows())
            << "First dimension differs: " << vec1.size() << " != " << vec2.rows();
        ASSERT_EQ(vec1[0].size(), vec2.cols())
            << "Second dimension differs: " << vec1[0].size() << " != " << vec2.cols();
    }

    TEST(NumpyArrayMethods, ArangeEigen){
        Eigen::VectorXd vec1(5);
        vec1 << 1., 2., 5., 9., 15.;
        // Vector after np.diff
        Eigen::VectorXd vec2(4);
        vec2 << 1., 3., 4., 6.;
        Eigen::VectorXd result = tsutil::adjDiffEigen(vec1);
        EXPECT_TRUE(compareVectors(vec2, result));
    }
}