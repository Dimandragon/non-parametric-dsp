#include <complex>
#include <cmath>
#include <vector>
#include <matplot/matplot.h>
#include <iostream>

import utility_math;


void createFill(std::vector<std::complex<float>> & data,
                 std::vector<std::pair<std::pair<float, float>, std::pair<float, float>>> f_series) {
    for (auto i = 0; i < data.size(); i++) {
        for (auto j = 0; j < f_series.size(); j++) {
            data[i] += f_series[j].first.first * std::cos(2.0 * std::numbers::pi / static_cast<double>(data.size()) * i * f_series[j].first.second);
            data[i] += f_series[j].second.first * std::sin(2.0 * std::numbers::pi / static_cast<double>(data.size()) * i * f_series[j].second.second);
        }
    }
}

void convertFSeriesFromComplexToTrigonometric(std::vector<std::complex<float>> compex_series, std::vector<std::pair<float, float>> & triginimetric_series) {
    triginimetric_series.clear();
    for (auto i = 0; i < compex_series.size(); i++) {
        triginimetric_series.push_back(
            std::pair{std::sqrt(compex_series[i].real()*compex_series[i].real() + compex_series[i].imag()*compex_series[i].imag()),
            std::atan(compex_series[i].imag()/compex_series[i].real())});
    }
}

struct Plotter {
    int plot_id = 0;

    void plotTrigonometricFourierSeries(std::vector<std::pair<float, float>> series, int len) {
        for (auto i = 0; i < series.size(); i++) {
            plot_id++;
            std::vector<float> plotting_data;
            for (auto j = 0; j < len; j++) {
                plotting_data.push_back(
                    std::round(series[i].first * std::cos( static_cast<float>(j) * static_cast<float>(i) * std::numbers::pi * 2.0  / static_cast<float>(len)
                    + series[i].second)*1000) / 1000);
            }

            matplot::plot(plotting_data);
            std::string filename = "/home/dmitry/projects/non-parametric-dsp/examples/utility/fft/images/triginometric_fourier_series";

            std::stringstream ss;
            ss << i;
            filename+=ss.str();
            filename+= ".svg";
            matplot::save(filename);

        }
        //atplot::save("/home/dmitry/projects/non-parametric-dsp/examples/utility/fft/images/triginometric_fourier_series.svg");
    }
};


int main(){
    auto len = 5000;
    std::vector<std::complex<float>> data(len);
    std::vector<std::complex<float>> out(len);

    std::vector<std::pair<std::pair<float, float>, std::pair<float, float>>> fill;
    fill.push_back({{1.0, 1.0}, {1.0, 1.0}});
    fill.push_back({{0.0, 0.0}, {5.0, 2.0}});
    fill.push_back({{10.0, 5.0}, {0.0, 0.0}});
    //fill.push_back({{5.0, 2.0}, {10.0, 12.0}});
    //createFill(data, fill);
    for (int i = 0; i < data.size(); i++) {
        //data[i] = static_cast<double>(i * 5) / static_cast<double>(data.size()) * std::cos(2.0 * std::numbers::pi / static_cast<double>(data.size()) * std::pow(static_cast<float>(i), 2.) / 10);
        data[i] = std::rand();
    }


    std::vector<float> plotting_data;
    for (int i = 0; i < data.size(); i++) {
        plotting_data.push_back(data[i].real());
    }
    //matplot::subplot(2, 2, 1);
    matplot::plot(plotting_data);
    matplot::save("/home/dmitry/projects/non-parametric-dsp/examples/utility/fft/images/data.svg");

    NP_DSP::ONE_D::UTILITY_MATH::fftc2c(data, out);

    std::vector<float> plotting_res_real;
    std::vector<float> plotting_res_imag;

    std::cout << out.size() << std::endl;

    for (auto i = 0; i < out.size(); i++){
        plotting_res_real.push_back(out[i].real());
        plotting_res_imag.push_back(out[i].imag());
    }

    matplot::plot(plotting_res_real);
    matplot::hold(matplot::on);
    matplot::plot(plotting_res_imag);
    matplot::hold(matplot::off);

    matplot::save("/home/dmitry/projects/non-parametric-dsp/examples/utility/fft/images/ft.svg");

    /*
    matplot::subplot(2, 2, 2);
    matplot::plot(plotting_res_real);
    matplot::hold(matplot::on);
    matplot::plot(plotting_res_imag);
    matplot::hold(matplot::off);
    */

    std::vector<std::pair<float, float>> trigonometric_series;
    convertFSeriesFromComplexToTrigonometric(out, trigonometric_series);
    std::vector<float> phase_vec;
    std::vector<float> ampl_vec;
    for (auto i = 0; i < trigonometric_series.size(); i++){
        phase_vec.push_back(trigonometric_series[i].second);
        ampl_vec.push_back(trigonometric_series[i].first);
    }
    matplot::plot(phase_vec);
    matplot::show();

    matplot::plot(ampl_vec);
    matplot::show();

    Plotter plotter;
    plotter.plot_id = 3;
    //plotter.plotTrigonometricFourierSeries(trigonometric_series, out.size());
    

    return 0;

}