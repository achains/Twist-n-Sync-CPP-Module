//
// Created by achains on 18.07.2021.
//

#include "TimeSync.h"

#include "xtensor/xnorm.hpp"

#include <algorithm>

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

    if (index_init > 0) {

    } else if (index_init < 0){

    } else {

    }
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