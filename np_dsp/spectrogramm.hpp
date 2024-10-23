#pragma once

/*
#include <cmath>
#include <matplot/matplot.h>
#include <set>
#include <thread>
#include <vector>
#include <functional>
#include <math.h>


void SpectrPlot(std::function<double(double, double)> data, size_t x_size, size_t y_size, double x_step, double y_step){
    using namespace matplot;

    //tiledlayout(2, 2);
    
    std::vector<double> x;
    std::vector<double> y;
    for (int i = 0; i < x_size; i++){
        x.push_back(i * x_step);
    }
    for (int i = 0; i < y_size; i++){
        y.push_back(i * y_step);
    }

    matplot::vector_2d Z(x.size(), matplot::vector_1d(y.size(), 0.));
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < y.size(); ++j) {
            auto x__ = x[i];
            auto y__ = y[i];
            Z[i][j] = data(x__, y__);
        }
    }

    imagesc(Z);

    show();
}

template<typename ExtractorT>
void spectrogramm(ExtractorT & extractor, size_t x_size, size_t y_size){
    std::vector<std::vector<double>> data;
    for (int i = 0; i < y_size; i++){
        data.push_back(std::vector<double>{});
        for (int j = 0; j < x_size; j++){
            data[i].push_back(0.0);
        }
    }
    auto freqToYIdx = [&](double freq) -> int{
        auto my_log = [&](double arg){
            return log2(arg) * log2(arg);
        };
        
        int idx = y_size - my_log(freq * extractor.getDataSize()) / my_log(0.5 * extractor.getDataSize()) * (y_size - 1) - 1;
        if (idx < 0){
            idx = 0;
        }
        if (idx > y_size - 1){
            idx = y_size - 1;
        }
        return idx;
    };

    auto idxToXIdx = [&](double idx) -> int {
        return idx / extractor.getDataSize() * (x_size - 1);
    };
    auto amplToZ = [&](double ampl) -> double{
        return std::log10(ampl);
    };

    for (int i = 0; i < extractor.getModesCount(); i++){
        std::vector<double> ampl_vec = extractor.getInstAmpl(i);
        std::vector<double> freq_vec = extractor.getInstFreq(i);
        for (int j = 0; j < extractor.getDataSize(); j++){
            data[freqToYIdx(freq_vec[j])][idxToXIdx(j)] = amplToZ(ampl_vec[j]);
        }
    }

    matplot::image(data, true);
    matplot::colorbar();

    matplot::show();
}

*/