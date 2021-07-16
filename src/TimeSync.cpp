//
// Created by achains on 18.07.2021.
//

#include "TimeSync.h"

#include "unsupported/Eigen/FFT"
#include "unsupported/Eigen/Splines"

#include <numeric>
#include <algorithm>
#include <memory>

TimeSync::TimeSync(Eigen::MatrixX3d  gyro_first, Eigen::MatrixX3d  gyro_second,
                   Eigen::VectorXd  ts_first, Eigen::VectorXd  ts_second, bool do_resample=true):
                        gyro_first_(std::move(gyro_first)), gyro_second_(std::move(gyro_second)),
                        ts_first_(std::move(ts_first)), ts_second_(std::move(ts_second)),
                        do_resample_(do_resample) {}


Eigen::ArrayXd TimeSync::adjDiffEigen(const Eigen::ArrayXd &A) {
    std::unique_ptr<Eigen::ArrayXd> B = std::make_unique<Eigen::ArrayXd>(A.size());
    std::adjacent_difference(A.begin(), A.end(), B->begin());

    return std::move(*B);
}

Eigen::VectorXd TimeSync::arangeEigen(const double &start, const double &finish, const double &step) {
    return Eigen::VectorXd::LinSpaced(static_cast<long>((finish - start) / step), start, finish);
}

void TimeSync::resample(double const & accuracy){
    double gyro_first_mean = adjDiffEigen(gyro_first_).mean(); // Тут надо взять [1: ], но потом ))
    double gyro_second_mean = adjDiffEigen(gyro_second_).mean();

    double dt = std::min({accuracy, gyro_first_mean, gyro_second_mean});

    if (do_resample_){
        Eigen::VectorXd ts_first_new = TimeSync::arangeEigen(ts_first_[0], *ts_first_.end(), dt);
        Eigen::VectorXd ts_second_new = TimeSync::arangeEigen(ts_second_[0], *ts_second_.end(), dt);

        //TODO: Interpolation function
    }
}

void TimeSync::obtainDelay(){
    // Correction of index numbering
    size_t shift = -gyro_first_.rows() + 1;
    // Cross-cor estimation
    std::vector<double> silence_vec_1;
    std::vector<double> silence_vec_2;

    std::vector<double> cross_cor = TimeSync::eigenCrossCor(silence_vec_1, silence_vec_2);
    size_t index_init = std::distance(cross_cor.begin(), std::max(cross_cor.begin(), cross_cor.end()));

    Eigen::MatrixX3d tmp_xx1;
    Eigen::MatrixX3d tmp_xx2;

    if (index_init > 0){
        1;
    }
    else if (index_init < 0){
        2;
    }

    size_t size = std::min(gyro_first_.rows(), gyro_second_.rows());
    tmp_xx1 = tmp_xx1(Eigen::seq(0, size - 1), Eigen::all);
    tmp_xx2 = tmp_xx2(Eigen::seq(0, size - 1), Eigen::all);

    // Calibration
    Eigen::MatrixX3d M = (tmp_xx1.transpose() * tmp_xx2) * (tmp_xx1.transpose() * tmp_xx1).inverse();

    // Cross-correlation re-estimation
    cross_cor = TimeSync::eigenCrossCor(silence_vec_1, silence_vec_2);
    index_init = std::distance(cross_cor.begin(), std::max(cross_cor.begin(), cross_cor.end()));

    // Cross-cor, based cubic spline coefficients
}


std::vector<double>
TimeSync::eigenCrossCor(std::vector<double> & xCorrInputVecFirst, std::vector<double> & xCorrInputVecSecond) {
    Eigen::FFT<double> fft;
    int N = std::max(xCorrInputVecFirst.size(), xCorrInputVecSecond.size());

    //Compute the FFT size as the "next power of 2" of the input vector's length (max)
    int b = ceil(log2(2.0 * N - 1));
    int fftsize = pow(2,b);
    int end = fftsize - 1;
    int maxlag = N - 1;
    size_t firstSize = xCorrInputVecFirst.size();
    size_t secondSize = xCorrInputVecSecond.size();

    //Zero Padd
    for (int i = xCorrInputVecFirst.size(); i < fftsize; i++)
    {
        xCorrInputVecFirst.push_back(0);
    }

    for (int i = xCorrInputVecSecond.size(); i < fftsize; i++)
    {
        xCorrInputVecSecond.push_back(0);
    }

    std::vector<std::complex<double> > freqvec;
    std::vector<std::complex<double> > freqvec2;

    //FFT for freq domain to both vectors
    fft.fwd( freqvec,xCorrInputVecFirst);
    fft.fwd( freqvec2,xCorrInputVecSecond);

    //Main step of cross corr
    for (int i = 0; i < fftsize; i++)
    {
        freqvec[i] = freqvec[i] * std::conj(freqvec2[i]);
    }

    std::vector<double> result;
    fft.inv(result, freqvec);

    //Will get rid of extra zero padding and move minus lags to beginning without copy
    std::vector<double> result2(std::make_move_iterator(result.begin() + end - maxlag + 1),
                                std::make_move_iterator(result.end()));

    result2.insert(result2.end(), make_move_iterator(result.begin())
            , make_move_iterator(result.begin()+maxlag));


    auto minMaxRange = std::minmax_element(result2.begin(),
                                           result2.end());

    //Will take back the changes which made in input vector
    if (xCorrInputVecFirst.size() != firstSize)
        xCorrInputVecFirst.resize(firstSize);
    if (xCorrInputVecSecond.size() != secondSize)
        xCorrInputVecSecond.resize(secondSize);

    return result2;
}

