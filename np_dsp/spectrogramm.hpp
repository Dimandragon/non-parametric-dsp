#pragma once

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
    //auto [X, Y] = meshgrid(x, y);

    //IC(x.size(), y.size());

    //const size_t n_rows = std::min(X.size(), Y.size());
    //const size_t n_cols = std::min(X[0].size(), Y[0].size());
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
        /*IC(log2(freq * extractor.getDataSize()), log2(0.5 * extractor.getDataSize()), 
            log2(freq * extractor.getDataSize()) / log2(0.5 * extractor.getDataSize()) * (y_size - 1));*/
        auto my_log = [&](double arg){
            return log2(arg) * log2(arg);
        };
        
        int idx = y_size - my_log(freq * extractor.getDataSize()) / my_log(0.5 * extractor.getDataSize()) * (y_size - 1) - 1;
        //int idx = (0.5 / freq / extractor.getDataSize() * 2.0 * y_size); 
        if (idx < 0){
            idx = 0;
            //std::cout << "small idx" << std::endl;
            //IC(my_log(freq * extractor.getDataSize()), my_log(0.5 * extractor.getDataSize()), 
            //    my_log(freq * extractor.getDataSize()) / my_log(0.5 * extractor.getDataSize()) * (y_size - 1));
        }
        if (idx > y_size - 1){
            idx = y_size - 1;
            //std::cout << "large idx" << std::endl;
            //IC(my_log(freq * extractor.getDataSize()), my_log(0.5 * extractor.getDataSize()), 
            //    my_log(freq * extractor.getDataSize()) / my_log(0.5 * extractor.getDataSize()) * (y_size - 1));
        }
        return idx;
        //return freq / 0.5 * (y_size - 1);
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

            //IC(j, idxToXIdx(j), freq_vec[j], freqToYIdx(freq_vec[j]), ampl_vec[j]);
            data[freqToYIdx(freq_vec[j])][idxToXIdx(j)] = amplToZ(ampl_vec[j]);
        }
    }

    /*std::function<double(double, double)> plotLamda = [&](double x, double y) -> double{
        return data[(int)x][(int)y];
    };

    SpectrPlot(plotLamda, x_size, y_size, (x_size - 1) / extractor.getDataSize(), (y_size - 1) / 0.5);*/
    matplot::image(data, true);
    matplot::colorbar();

    matplot::show();
}