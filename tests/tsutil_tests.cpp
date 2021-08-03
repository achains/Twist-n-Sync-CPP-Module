//
// Created by achains on 03.08.2021.
//


#include "util/TSUtil.h"
#include "tsutil_tests.h"

namespace tsutil_tests{
    template <typename T1, typename T2>
    testing::AssertionResult compareVectors(T1 const & vec1, T2 const & vec2, double eps){

        if (vec1.size() != vec2.size())
            return testing::AssertionFailure() << "vec1.size() != vec2.size() <<>> "
                                               << vec1.size() << " != " << vec2.size();

        for (Eigen::Index i = 0; i < vec1.size(); ++i)
            if (std::abs(vec1[i] - vec2[i]) > eps)
                return testing::AssertionFailure() << "Elements at index " << i << " are not equal with eps=" << eps;
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

    // np.diff()
    TEST(NumpyArrayMethods, AdjDiffEigen){
        Eigen::VectorXd vec1(5);
        vec1 << 1., 2., 5., 9., 15.;
        // Vector after np.diff
        Eigen::VectorXd expected(4);
        expected << 1., 3., 4., 6.;
        Eigen::VectorXd result = tsutil::adjDiffEigen(vec1);
        EXPECT_TRUE(compareVectors(expected, result));
    }

    // np.arange()
    TEST(NumpyArrayMethods, ArangeEigen){
        // np.arange(1.0, 10.0, 2.0)
        Eigen::VectorXd vec1(5);
        vec1 << 1., 3., 5., 7., 9.;
        Eigen::VectorXd vec1_result = tsutil::arangeEigen(1.0, 10.0, 2.0);
        EXPECT_TRUE(compareVectors(vec1, vec1_result));

        // np.arange(1.0, 5.0)
        Eigen::VectorXd vec2(4);
        vec2 << 1., 2., 3., 4;
        Eigen::VectorXd vec2_result = tsutil::arangeEigen(1.0, 5.0);
        EXPECT_TRUE(compareVectors(vec2, vec2_result));
    }

    // np.linalg.norm(matrix : np.array(N, 3), axis=1)
    TEST(NumpyArrayMethods, NormOfRows){
        Eigen::MatrixX3d vec1(4, 3);
        vec1 << 3., 16., 5.,
                5., 6., 11.,
                13., 12., 10.,
                1., 2., 3.;
        // Norms calculated by np.linalg.norm
        Eigen::VectorXd expected(4);
        expected << 17.02938637, 13.49073756, 20.32240143, 3.74165739;
        Eigen::VectorXd result = tsutil::getNormOfRows(vec1);
        EXPECT_TRUE(compareVectors(expected, result, 1e-8));
    }

    // np.roots()
    TEST(QuadraticEquation, QuadraticSolver){
        Eigen::VectorXd coefficients(3);
        coefficients << 1., 4., 4.;
        Eigen::VectorXd expected(2);
        expected << -2., -2.;
        Eigen::VectorXd result = tsutil::quadraticRoots(coefficients.reverse());
        EXPECT_TRUE(compareVectors(expected, result, 1e-8));
    }

    // scipy.signal.correlate
    TEST(CrossCorrelation, CrossCorrEqualSize){
        Eigen::VectorXd vec1(5);
        Eigen::VectorXd vec2(5);
        vec1 << 3., 17., 19., 12., 6.;
        vec2 << 4., 3., 1., 13., 7.;
        Eigen::VectorXd expected(9);
        expected << 21., 158., 357., 357., 280., 215., 118.,  66.,  24.;
        Eigen::VectorXd result = tsutil::eigenCrossCor(vec1, vec2);
        EXPECT_TRUE(compareVectors(expected, result, 1e-8));
    }

    // scipy.signal.correlate
    TEST(CrossCorrelation, CrossCorrRandomSize){
        Eigen::VectorXd vec1(3);
        Eigen::VectorXd vec2(5);
        vec1 << 3., 1.16, 17.3;
        vec2 << 4., 5.6, 11.2, 33.2, 1.5;
        Eigen::VectorXd expected(7);
        expected << 4.5, 101.34, 98.062, 604.152, 212.256, 101.52, 69.2;
        Eigen::VectorXd result = tsutil::eigenCrossCor(vec1, vec2);
        EXPECT_TRUE(compareVectors(expected, result, 1e-8));
    }
}