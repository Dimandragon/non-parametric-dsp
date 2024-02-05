//
// Created by dmitry on 21.12.23.
//
#include <cmath>
#include <vector>
#include <matplot/matplot.h>

void fill(std::vector<double> & data, double ampl) {
    for (auto i = 0; i < data.size(); i++) {
        data[i] = ampl * std::cos(2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(data.size()));
    }
}

double loss(std::vector<double> data1, std::vector<double> data2) {
    double accum = 0.0;
    for (auto i = 0; i < data1.size(); i++){
        accum += (data1[i] - data2[i]) * (data1[i] - data2[i]);
    }
    return accum;
}

int main(){
    auto len = 128;
    std::vector<double> data1(len);
    std::vector<double> data2(len);
    std::vector<double> losses;

    fill(data1, 1.0);

    for (auto i = 0; i < data1.size(); i++) {
        data1[i] = data1[i];
    }
    matplot::plot(data1);
    matplot::save("/home/dmitry/projects/non-parametric-dsp/examples/utility/fft/images/cos_first.svg");

    for (auto i = 0; i < data1.size(); i++) {
        fill(data2, static_cast<double>(i) * 2.0 / static_cast<double>(data1.size()));
        matplot::plot(data2);
        std::string filename = "/home/dmitry/projects/non-parametric-dsp/examples/utility/fft/images/cos_second";

        std::stringstream ss;
        ss << i;
        filename+=ss.str();
        filename+= ".svg";
        matplot::save(filename);

        losses.push_back(loss(data1, data2));
    }
    matplot::plot(losses);
    matplot::save("/home/dmitry/projects/non-parametric-dsp/examples/utility/fft/images/losses.svg");
    return 0;
}