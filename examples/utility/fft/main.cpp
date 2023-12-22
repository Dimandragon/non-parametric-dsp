#include <complex>
#include <cmath>
#include <vector>
#include "pocketfft_hdronly.h"
#include <matplot/matplot.h>

import utility_math;
import <iostream>;

//using namespace std;
using namespace pocketfft;

// floating point RNG which is good enough for simple demos
// Do not use for anything important!
inline double simple_drand()
{
    constexpr double norm = 1./RAND_MAX;
    return rand()*norm;
}

template<typename T> void crand(std::vector<std::complex<T>> &v)
{
    for (auto & i:v)
        i = std::complex<T>(simple_drand()-0.5, 0.);
}

template<typename T> void vec_sin(std::vector<std::complex<T>> &v)
{
    for (int i = 0; i < v.size(); i++)
        v[i] = std::complex<T>(sinf(((float)i)) - 1.0, 0.f);
}

template<typename T> void crand(std::vector<T> &v)
{
    for (auto & i:v)
        i = simple_drand()-0.5;
}

template<typename T1, typename T2> long double l2err
        (const std::vector<T1> &v1, const std::vector<T2> &v2)
{
    long double sum1=0, sum2=0;
    for (size_t i=0; i<v1.size(); ++i)
    {
        long double dr = v1[i].real()-v2[i].real(),
                di = v1[i].imag()-v2[i].imag();
        long double t1 = std::sqrt(dr*dr+di*di), t2 = std::abs(v1[i]);
        sum1 += t1*t1;
        sum2 += t2*t2;
    }
    return std::sqrt(sum1/sum2);
}


int main()
{
    size_t len = 100;

    shape_t shape{len};
    stride_t stridef(shape.size());
    size_t tmpf=sizeof(std::complex<float>);
    for (int i=shape.size()-1; i>=0; --i)
    {
        stridef[i]=tmpf;
        tmpf*=shape[i];
    }
    size_t ndata=1;
    for (size_t i=0; i<shape.size(); ++i)
    {
        ndata*=shape[i];
    }

    std::vector<std::complex<float>> dataf(len);
    std::vector<std::complex<float>> dataf2(len);
    vec_sin(dataf);
    vec_sin(dataf2);

    shape_t axes;
    for (size_t i=0; i<shape.size(); ++i)
    {
        axes.push_back(i);
    }

    auto resf = dataf;
    auto resf2 = dataf2;

    std::vector<float> plotting_vector_first_real;
    std::vector<float> plotting_vector_first_imag;
    for (auto i = 0; i < 100; i++){
        plotting_vector_first_real.push_back(dataf[i].real());
        plotting_vector_first_imag.push_back(dataf[i].imag());
    }

    matplot::subplot(2, 2, 0);
    matplot::plot(plotting_vector_first_real);
    matplot::hold(matplot::on);
    matplot::plot(plotting_vector_first_imag);

    c2c(shape, stridef, stridef, axes, FORWARD,
        dataf.data(), resf.data(), (float)0.01);
    NP_DSP::ONE_D::UTILITY_MATH::fftc2c(dataf2, resf2);

    for (auto i = 0; i < resf.size(); i++){
        //std::cout << dataf[i] << " " << dataf2[i] << "    ";
    }
    std::cout << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << std::endl;
    for (auto i = 0; i < resf.size(); i++){
        std::cout << resf[i] << " " << resf2[i] << "    ";
    }
    std::cout << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << std::endl;

    std::vector<float> plotting_vector_real{};
    std::vector<float> plotting_vector_imag{};
    for (auto i = 0; i < 100; i++){
        plotting_vector_real.push_back(0.);
        plotting_vector_imag.push_back(0.);
    }

    for (auto i = 0; i < resf.size(); i++){
        plotting_vector_real[i] = resf[i].real();
        plotting_vector_imag[i] = resf[i].imag();
    }
    matplot::subplot(2, 2, 1);
    matplot::plot(plotting_vector_real);
    matplot::hold(matplot::on);
    matplot::plot(plotting_vector_imag);
    matplot::hold(matplot::off);

    c2c(shape, stridef, stridef, axes, BACKWARD,
        resf.data(), dataf.data(), 1.f);
    NP_DSP::ONE_D::UTILITY_MATH::fftc2c(resf2, dataf2);

    for (auto i = 0; i < resf.size(); i++){
        std::cout << resf[i] << " " << resf2[i] << "    ";
    }
    std::cout << std::endl;


    std::vector<float> plotting_vector_last_real;
    std::vector<float> plotting_vector_last_imag;

    for (auto i = 0; i < 100; i++){
        plotting_vector_last_real.push_back(dataf[i].real());
        plotting_vector_last_imag.push_back(dataf[i].imag());
    }
    matplot::subplot(2, 2, 2);
    matplot::plot(plotting_vector_last_real);
    matplot::hold(matplot::on);
    matplot::plot(plotting_vector_last_imag);
    matplot::hold(matplot::off);
    matplot::save("/home/dmitry/projects/non-parametric-dsp/examples/utility/fft/images/data.svg");
    return 0;
}
