#include <iostream>
#include <algorithm>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xadapt.hpp>

#include <xtensor-blas/xlinalg.hpp>


using namespace xt::placeholders;

int main() {
    xt::xarray<double> test_vec = {{1., 2., 3., 4.}, {4., 5., 6., 7.}, {8., 2., 4., 7.}};
    xt::xarray<double> next_vec = {{0., 0., 1.}, {1., 0., 0.}, {0., 1., 0.}};
    xt::xarray<double> vec_two{{1., 2.}};

    std::cout << xt::norm_l2(test_vec, {0}) << std::endl;

    std::cout << test_vec.shape()[1] << std::endl;

    std::cout << xt::linalg::inv(next_vec);

    return 0;
}