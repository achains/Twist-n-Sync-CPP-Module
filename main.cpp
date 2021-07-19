#include <iostream>
#include <algorithm>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xadapt.hpp>

int main() {
    xt::xarray<double> test_vec = {{1., 2., 3.}, {4., 5., 6.}};
    xt::xarray<double> vec_two{{1., 2.}};

    std::cout << xt::norm_l2(test_vec, {0}) << std::endl;

    std::cout << test_vec.shape()[1];

    return 0;
}