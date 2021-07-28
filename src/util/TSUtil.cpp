//
// Created by achains on 23.07.2021.
//

#include "TSUtil.h"
#include "unsupported/Eigen/FFT"
#include "util/CubicSpline.h"

#include <numeric>

namespace TSUtil {

    Eigen::VectorXd arangeEigenDeprecate(const double &start, const double &finish, const double &step) {
        return Eigen::VectorXd::LinSpaced(static_cast<long>(std::ceil((finish - start) / step)), start, finish);
    }

    Eigen::VectorXd arangeEigen(double start, double const & stop, double const & step){
        Eigen::Index size = std::ceil((stop - start) / step);
        Eigen::VectorXd result(size);
        for (Eigen::Index i = 0; i < size; ++i) {
            result[i] = start;
            start += step;
        }
        return result;
    }

    Eigen::VectorXd adjDiffEigen(const Eigen::VectorXd &data) {
        Eigen::VectorXd new_data(data.size());
        std::adjacent_difference(data.begin(), data.end(), new_data.begin());

        return new_data(Eigen::seq(1, Eigen::last));
    }

    Eigen::VectorXd vectorToEigVectorXd(std::vector<double> &data) {
        return Eigen::Map<Eigen::VectorXd>(data.data(), static_cast<int> (data.size()));
    }

    Eigen::MatrixX3d vectorToEigMatrixX3d(std::vector<std::vector<double>> &data) {
        Eigen::MatrixX3d eigen_data(data.size(), 3);
        for (int i = 0; i < data.size(); ++i)
            eigen_data.row(i) = Eigen::Map<Eigen::VectorXd>(data[i].data(), 3);
        return eigen_data;
    }

    Eigen::VectorXd getNormOfRows(Eigen::MatrixX3d const & data){
        Eigen::VectorXd norm_data(data.rows());
        Eigen::Index i = 0;
        for (auto row : data.rowwise()){
            norm_data[i++] = row.norm();
        }
        return norm_data;
    }

    Eigen::VectorXd interpolate(Eigen::VectorXd const & x_old, Eigen::VectorXd const & y_old,
                                Eigen::VectorXd const & x_new){

        CubicSpline interp(x_old, y_old);
        return interp.getValuesOnSegment(x_new);
    }

    std::vector<double> eigenCrossCor(std::vector<double> & data_1, std::vector<double> & data_2) {
        // Cross-cor(x, y) = iFFT(FFT(x) * conj(FFT(y)))

        long shift_size = static_cast<long>(data_1.size());

        // Length of Discrete Fourier Transform
        // TODO: It may be faster if we round N up to the next power of two (std::vector eigenCrossCor)
        size_t N = data_1.size() + data_2.size() - 1;

        // Zero padding both vectors
        data_1.resize(N);
        data_2.resize(N);

        Eigen::FFT<double> fft;

        std::vector<std::complex<double>> fft_first;
        std::vector<std::complex<double>> fft_second;

        fft.fwd(fft_first, data_1);
        fft.fwd(fft_second, data_2);

        std::vector<std::complex<double>> fft_result(N);
        for (size_t i = 0; i < N; ++i)
            fft_result[i] = fft_first[i] * std::conj(fft_second[i]);

        std::vector<double> cross_cor(N);
        fft.inv(cross_cor, fft_result);

        // Rotating cross-correlation vector on shift_size
        std::rotate(cross_cor.begin(), cross_cor.begin() + shift_size, cross_cor.end());

        return cross_cor;
    }

    Eigen::VectorXd eigenCrossCor(Eigen::VectorXd & data_1, Eigen::VectorXd & data_2){
        // Cross-cor(x, y) = iFFT(FFT(x) * conj(FFT(y)))

        Eigen::Index shift_size = data_1.size();

        // Length of Discrete Fourier Transform
        // TODO: It may be faster if we round N up to the next power of two (Eigen::VectorXd eigenCrossCor)
        Eigen::Index N = data_1.size() + data_2.size() - 1;

        // Zero padding both vectors
        data_1.conservativeResize(N);
        data_2.conservativeResize(N);

        Eigen::FFT <double> fft;

        Eigen::VectorXcd fft_first(N);
        Eigen::VectorXcd fft_second(N);

        fft.fwd(fft_first, data_1);
        fft.fwd(fft_second, data_2);


        Eigen::VectorXcd fft_result = fft_first.array() * fft_second.conjugate().array();

        Eigen::VectorXd unbiased_result(N);
        fft.inv(unbiased_result, fft_result);

        // Rotating cross-correlation vector on shift_size
        Eigen::VectorXd cross_cor(N);
        cross_cor << unbiased_result(Eigen::seq(shift_size, Eigen::last)),
                     unbiased_result(Eigen::seq(0, shift_size - 1));

        return cross_cor;
    }
}