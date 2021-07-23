#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/Splines>
#include <deque>
#include <vector>
#include <algorithm>
#include <libInterpolate/Interpolate.hpp>

#include "util/TSUtil.h"

int main()
{
    Eigen::MatrixX3d gyro_first(4, 3);
    gyro_first << 1., 1., 1.,
                  2., 2., 2.,
                  3., 3., 3.,
                  4., 4., 4.;


    Eigen::MatrixX3d gyro_second(6, 3);
    gyro_second << 2., 5., 3.,
                   4., 7., 8.,
                   1., 2., 2.5,
                   12.3, 11.5, 10.5,
                   5.3, 2., 1.,
                   8.2, 3.4, 1.2;

    Eigen::VectorXd norm_first = TSUtil::getNormOfRows(gyro_first);
    Eigen::VectorXd norm_second = TSUtil::getNormOfRows(gyro_second);

    std::cout << TSUtil::eigenCrossCor(norm_first, norm_second) << std::endl;

    // Equal to np [:-1]
    std::cout << gyro_first(Eigen::seq(0, Eigen::last - 1), Eigen::all);
    return 0;
}