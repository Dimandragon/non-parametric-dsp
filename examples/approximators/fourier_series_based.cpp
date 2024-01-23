#include <icecream.hpp>

import signals;
import approximators;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;
import <cstdlib>;
import <complex>;
import <iostream>;
import <matplot/matplot.h>;
import <string>;
import utility_math;

void createFill(auto & signal){
    IC(signal.size());
    for (auto i = 0; i < signal.size(); i++){
        IC(i);
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
    IC(incr(point)); //0
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>> signal1;
    using SignalT = decltype(signal1);
    IC(incr(point)); //1
    
    for (auto i = 0; i < 100; i++){
        signal1.base->vec->push_back(0.);
        IC(i);
    }
    createFill(signal1);
    IC(incr(point)); //2

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/approximators/images/signal1.svg");
    SignalT signal2;
    auto error = [&](auto & approximator){
        approximator.is_actual = false;
        auto accum = 0.;
        for(auto i = 0; i < signal1.size(); i++){
            accum += (approximator.compute(i) - signal1[i]) * (approximator.compute(i) - signal1[i]);
            //accum += std::abs(approximator.compute(i) - signal1[i]);
            //IC(std::abs(approximator.compute(i) - signal1[i]));
        }
        //approximator.show(NP_DSP::ONE_D::PlottingKind::Simple);
        return accum;
    };

    auto stopPoint = [](auto losses_different, auto & approximator) {
        if (losses_different > 0.01){ //todo move precision to external parameter
            return false;
        }
        else{
            return true;
        }
    };
    
    NP_DSP::ONE_D::APPROX::FourierSeriesBased<SignalT, decltype(error), decltype(stopPoint)> approximator(error, signal1, stopPoint);
    approximator.train();
    approximator.is_actual = false;
    for (auto i = 0; i < signal1.size(); i++){
        signal2.base->vec->push_back(approximator.compute(i));
    }
    
    approximator.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/approximators/images/signal2.svg");

    plotComplex(approximator.fourier_series, "/home/dmitry/projects/non-parametric-dsp/examples/approximators/images/approximator_fs.svg");
    std::vector<std::complex<double>> data_vec;
    std::vector<std::complex<double>> data_fs;
    for (auto i = 0; i < signal1.size(); i++){
        data_vec.push_back({signal1[i], 0.});
        data_fs.push_back({0., 0.});
    }
    NP_DSP::ONE_D::UTILITY_MATH::fftc2c(data_vec, data_fs);
    plotComplex(data_fs, "/home/dmitry/projects/non-parametric-dsp/examples/approximators/images/data_fs.svg");
    return 0;
}