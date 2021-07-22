#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/Splines>
#include <deque>
#include <vector>
#include <algorithm>

Eigen::VectorXf arange(const float& start, const float& finish, const float& step){
    return Eigen::VectorXf::LinSpaced(static_cast<long>((finish - start) / step), start, finish);
}

std::vector<double> eigenCrossCorrelation(std::vector<double>& xCorrInputVecFirst
        , std::vector<double>& xCorrInputVecSecond);

int main()
{
    std::vector<double> a = {1., 2., 3.};
    std::cout << arange(0., 5., 0.5) << std::endl;

    Eigen::VectorXf tmp_vec {{1., 2., 3.}};
    std::cout << *tmp_vec.begin() << ' ' << tmp_vec.tail(1) << std::endl;

    std::vector<double> time1 {1.1, 2.2, 3., 121.5};
    std::vector<double> time2 {1.5, 1.6, 1.5, 9.2, 5.8};
    std::vector<double> result = eigenCrossCorrelation(time1, time2);

    for (auto &elem: result){
        std::cout << elem << ' ';
    }

    std::cout << "\n\n\n\n";
    std::vector<double> std_from_eigen_vec(tmp_vec.begin(), tmp_vec.end());
    for (auto & elem: std_from_eigen_vec) std::cout << elem << ' ';

    std::cout << std::endl;
    Eigen::MatrixX3d test_dim_matrix {{1., 2., 3.}, {4., 5., 6.}, {4., 5., 6.}, {3., 2., 1.}};

    std::cout << test_dim_matrix.rows() << std::endl << std::endl;
    std::cout << "Slice Experiments: \n\n";
    std::cout << test_dim_matrix(Eigen::seq(0, 2), Eigen::all) << std::endl;


    // Experiments with Splines
    std::cout << "Spline Experiments: \n\n";

    Eigen::RowVectorXd time{{1., 2., 3., 4., 5.,}};
    Eigen::RowVectorXd signal{{5.7, 3.4, 5.7, 6.8, 10.5}};

    typedef Eigen::Spline<double, 1, 3> CubicSpline1D;
    typedef Eigen::SplineFitting<CubicSpline1D> CubicSplineFitting1D;

    // Generating spline function
    const auto fit = CubicSplineFitting1D::Interpolate(signal, 3, time);
    const CubicSpline1D& spline(fit);

    std::cout << "spline.ctrls(): " << spline.ctrls() << std::endl;
    std::cout << "spline.knots(): " << spline.knots() << std::endl;
    for (int i = 0; i < 10; ++i){
        std::cout << '(' << i << ", " << spline(i).coeff(0) << ')' << std::endl;
    }

    return 0;
}


std::vector<double> eigenCrossCorrelation(std::vector<double>& xCorrInputVecFirst
        , std::vector<double>& xCorrInputVecSecond)
{
    Eigen::FFT<double> fft;
    size_t N = std::max(xCorrInputVecFirst.size(), xCorrInputVecSecond.size());

    //Compute the FFT size as the "next power of 2" of the input vector's length (max)
    int b = ceil(log2(2.0 * static_cast<double>(N) - 1));
    auto fftsize = static_cast<size_t>(pow(2,b));
    size_t end = fftsize - 1;

    size_t maxlag = N - 1;
    size_t firstSize = xCorrInputVecFirst.size();
    size_t secondSize = xCorrInputVecSecond.size();

    //Zero Padd
    for (size_t i = xCorrInputVecFirst.size(); i < fftsize; ++i){
        xCorrInputVecFirst.push_back(0);
    }

    for (size_t i = xCorrInputVecSecond.size(); i < fftsize; ++i){
        xCorrInputVecSecond.push_back(0);
    }

    std::vector<std::complex<double> > freqvec;
    std::vector<std::complex<double> > freqvec2;

    //FFT for freq domain to both vectors
    fft.fwd( freqvec,xCorrInputVecFirst);
    fft.fwd( freqvec2,xCorrInputVecSecond);

    //Main step of cross corr
    for (size_t i = 0; i < fftsize; ++i)
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


    //Return val
    auto resultIndex = ((minMaxRange.second - result2.begin()) - N + 1);
    std::cout << "Max element at: "  << resultIndex << std::endl;
    auto maxValue = result[minMaxRange.second - result.begin()];
    return result2;
}