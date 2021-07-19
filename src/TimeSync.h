//
// Created by achains on 18.07.2021.
//

#ifndef TWIST_N_SYNC_CPP_MODULE_TIMESYNC_H
#define TWIST_N_SYNC_CPP_MODULE_TIMESYNC_H

#include "alglib-3.17.0/src/interpolation.h"
#include "alglib-3.17.0/src/fasttransforms.h"

#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"

#include <tuple>


class TimeSync {
 public:
    explicit TimeSync(xt::xarray<double> gyro_first, xt::xarray<double> gyro_second,
                      xt::xarray<double> ts_first, xt::xarray<double> ts_second, bool do_resample);
 private:
    // 3D angular velocities from devices' gyros
    xt::xarray<double> gyro_first_;
    xt::xarray<double> gyro_second_;

    // Gyros' timestamps
    xt::xarray<double> ts_first_;
    xt::xarray<double> ts_second_;

    // Resampled timestamps
    xt::xarray<double> ts_first_re_;
    xt::xarray<double> ts_second_re_;

    double dt_{};

    // Flag to do resampling of angular velocities
    bool do_resample_{};

    void resample(double const & accuracy);

    void obtain_delay();

    // Obtain initial index of argmax of cross-cor. Related to initial estimation of time delay
    xt::xarray<double> get_initial_index();

    xt::xarray<double> interpolate();
};


#endif //TWIST_N_SYNC_CPP_MODULE_TIMESYNC_H
