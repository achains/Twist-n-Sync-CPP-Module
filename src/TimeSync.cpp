//
// Created by achains on 18.07.2021.
//

#include "TimeSync.h"

#include "xtensor/xnorm.hpp"
#include "xtensor/xview.hpp"

#include <algorithm>

// Placeholder '_' in xt::range
using namespace xt::placeholders;

TimeSync::TimeSync(xt::xarray<double> gyro_first, xt::xarray<double> gyro_second,
                   xt::xarray<double> ts_first, xt::xarray<double> ts_second, bool do_resample=true):
                        gyro_first_(std::move(gyro_first)), gyro_second_(std::move(gyro_second)),
                        ts_first_(std::move(ts_first)), ts_second_(std::move(ts_second)), do_resample_(do_resample) {}

void TimeSync::resample(double const & accuracy) {
    dt_ = std::min({accuracy, xt::mean(xt::diff(ts_first_))(), xt::mean(xt::diff(ts_second_))()});

    if (do_resample_){
        ts_first_re_ = xt::arange<double>(ts_first_.front(), ts_first_.back() + dt_, dt_);
        ts_second_re_ = xt::arange<double>(ts_second_.front(), ts_second_.back() + dt_, dt_);

        // TODO: Finish interpolation function
        gyro_first_ = interpolate();
        gyro_second_ = interpolate();
    }
    else{
        ts_first_re_ = ts_first_;
        ts_second_re_ = ts_second_;
    }
}

xt::xarray<double> TimeSync::interpolate() {
    // Correction of index numbering
    size_t shift = -gyro_first_.shape()[0] + 1;
    // Cross-correlation estimation
    xt::xarray<double> cross_cor = get_initial_index();
    size_t index_init = std::distance(cross_cor.begin(), std::max_element(cross_cor.begin(), cross_cor.end())) + shift;

    xt::xarray<double> gyro_first_tmp;
    xt::xarray<double> gyro_second_tmp;

    if (index_init > 0) {
        gyro_first_tmp = xt::view(gyro_first_, xt::range(_, -index_init));
        gyro_second_tmp = xt::view(gyro_second_, xt::range(index_init, _));
    } else if (index_init < 0){
        gyro_first_tmp = xt::view(gyro_first_, xt::range(-index_init, _));
        gyro_second_tmp = xt::view(gyro_second_, xt::range(_, index_init));
    } else {
        gyro_first_tmp = gyro_first_;
        gyro_second_tmp = gyro_second_;
    }

    size_t right_border = std::min(gyro_first_tmp.shape()[0], gyro_first_tmp.shape()[1]);

    gyro_first_tmp = xt::view(gyro_first_tmp, xt::range(_, right_border));
    gyro_second_tmp = xt::view(gyro_second_tmp, xt::range(_, right_border));

    // TODO: Can't build xtensor-blas, undefined reference to `dgetrf_'
    // Calibration
    //xt::xarray<double> M = (xt::transpose(gyro_second_tmp) * gyro_first_tmp);

    // Cross-correlation re-estimation
    cross_cor = get_initial_index();
    index_init = std::distance(cross_cor.begin(), std::max_element(cross_cor.begin(), cross_cor.end()));

    // Cross-correlation, based cubic spline Coefficient
}

void TimeSync::obtain_delay() {

}

xt::xarray<double> TimeSync::get_initial_index() {
    xt::xarray<double> x1_temp = xt::norm_l2(gyro_first_, {1});
    xt::xarray<double> x2_temp = xt::norm_l2(gyro_second_, {1});
    // TODO: Fit cross-correlation function from alglib
    xt::xarray<double> cross_cor = {1.}; //alglib::corrr1d()

    return cross_cor;
}