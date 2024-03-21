#include <icecream.hpp>

#include <signals.hpp>
#include <approximators.hpp>
#include <cmath>
#include <vector>
#include <npdsp_concepts.hpp>
#include <string>
#include <cstdlib>
#include <complex>
#include <iostream>
#include <matplot/matplot.h>
#include <string>
#include <utility_math.hpp>

void createFill(auto & signal){
    for (auto i = 0; i < signal.size(); i++){
        signal[i] = std::rand();
    }
}

void plotComplex(auto & signal, const std::string & filename){
    std::vector<double> real;
    std::vector<double> imag;
    for (auto i = 0; i < signal.size(); i++){
        real.push_back(signal[i].real());
        imag.push_back(signal[i].imag());
    }
    matplot::plot(real);
    matplot::hold(matplot::on);
    matplot::plot(imag, "--");
    matplot::hold(matplot::off);
    matplot::show();
    matplot::save(filename);
}

int incr(int & point){
    point++;
    return point - 1;
}

int main(){
    int point = 0;
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    //if (len == 100) fs[51].re == fs[49].re, fs[51].im = -fs[49].im; fs[50]is uniq
    //if len == 101 its 50 and 51
    //0 and last are different
    //if 10 then 4 and 6
    //if 11 then 5 and 6
    for (auto i = 0; i < 102; i++){
        signal1.base->vec->push_back(0.);
    }
    createFill(signal1);

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/approximators/images/signal1.svg");
    SignalT signal2;
    auto error = [&](auto & approximator){
        approximator.is_actual = false;
        auto accum = 0.0;
        for(auto i = 0; i < signal1.size(); i++){
            accum += (approximator.template compute<double, size_t>(i) - signal1[i]) *
                    (approximator.template compute<double, size_t>(i) - signal1[i])/signal1.size();
        }
        return std::abs(accum);
    };

    auto bySampleError = [&](auto & approximator, auto i){
        return (approximator.approximated_data[i].real() - signal1[i]) * (approximator.approximated_data[i].real() - signal1[i])/ approximator.tile_size;
    };

    auto stopPoint = [](auto losses_different, auto & approximator) {
        if (losses_different > 0.0001){ //todo move precision to external parameter
            return false;
        }
        else{
            return true;
        }
    };

    NP_DSP::ONE_D::APPROX::FourierSeriesBased<decltype(error), decltype(stopPoint),
        NP_DSP::ONE_D::APPROX::FSApproxKind::Simple, decltype(bySampleError)> approximator(error, signal1, stopPoint);
    approximator.setApproxOrderRatio(1.);
    approximator.tile_size = 5;
    approximator.bySampleLoss = &bySampleError;
    approximator.train();
    approximator.is_actual = false;
    for (auto i = 0; i < signal1.size(); i++){
        signal2.base->vec->push_back(approximator.compute<double, size_t>(i));
    }

    
    approximator.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/approximators/images/signal2.svg");

    plotComplex(approximator.fourier_series, "/home/dmitry/projects/non-parametric-dsp/examples/approximators/images/approximator_fs.svg");
    std::vector<std::complex<double>> data_vec;
    std::vector<std::complex<double>> data_fs;
    
    double average = 0;
    for (auto i = 0; i < signal1.size(); i++){
        average += signal1[i] / signal1.size();
        data_vec.push_back({signal1[i], 0.});
        data_fs.push_back({0., 0.});
    }
    NP_DSP::ONE_D::UTILITY_MATH::fftc2c(data_vec, data_fs);
    plotComplex(data_fs, "/home/dmitry/projects/non-parametric-dsp/examples/approximators/images/data_fs.svg");
    return 0;
}