#include <iostream>
#include <matplot/matplot.h>
#include <set>
#include <thread>
#include <vector>
#include <icecream.hpp>

int main() {
    using namespace matplot;

    //tiledlayout(2, 2);
    size_t N = 200;
    
    std::vector<double> x = linspace(-1., +1., N);
    std::vector<double> y = linspace(-1., +1., N);
    auto [X, Y] = meshgrid(x, y);

    IC(x.size(), y.size());

    const size_t n_rows = std::min(X.size(), Y.size());
    const size_t n_cols = std::min(X[0].size(), Y[0].size());
    matplot::vector_2d Z(n_rows, matplot::vector_1d(n_cols, 0.));
    for (size_t i = 0; i < n_rows; ++i) {
        for (size_t j = 0; j < n_cols; ++j) {
            Z[i][j] = peaks(X[i][j], Y[i][j]);
        }
    }

    imagesc(Z);

    show();
    return 0;
}