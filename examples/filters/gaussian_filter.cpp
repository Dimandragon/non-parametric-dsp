#include <filters.hpp>
#include <signals.hpp>
#include <utility_math.hpp>

using namespace NP_DSP::ONE_D;

int main(){
    GenericSignal<SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT mask;
    std::vector<std::complex<double>> spectre1;
    std::vector<std::complex<double>> spectre2;
    std::vector<std::complex<double>> buffer1;
    std::vector<std::complex<double>> buffer2;
    std::vector<std::complex<double>> afr;

    size_t size = 500;

    for (int i = 0; i < size; i++){
        //mask.base->vec->push_back(0.0);
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0.0);
    }

    FILTERS::generateConvMask(0.1, 1.0, mask, FILTERS::MaskKind::Gaussian, 1);

    mask.show(PlottingKind::Simple);

    FILTERS::MonoFreqFilters<double,
        NP_DSP::ONE_D::FILTERS::MonoInstFreqFilteringType::Conv> filter;

    filter.conv_filter_len = mask.size();

    //filter.show(PlottingKind::Spectre);

    filter.compute(signal1, signal2, mask);

    signal1.show(PlottingKind::Simple);
    signal2.show(PlottingKind::Simple);

    for(int i = 0; i < signal1.size(); i++){
        buffer1.push_back({signal1[i], 0.0});
        buffer2.push_back({signal2[i], 0.0});
        spectre1.push_back({0.0, 0.0});
        spectre2.push_back({0.0, 0.0});
    }

    UTILITY_MATH::fftc2c(buffer1, spectre1);
    UTILITY_MATH::fftc2c(buffer2, spectre2);

    for(int i = 0; i < signal1.size(); i++){
        afr.push_back(spectre2[i] / spectre1[i]);
    }
    std::vector<double> plotting_vector;
    for (int i = 0; i < signal1.size(); i++){
        plotting_vector.push_back(afr[i].real());
    }
    matplot::plot(plotting_vector);
    matplot::hold(true);
    for (int i = 0; i < signal1.size(); i++){
        plotting_vector[i] = afr[i].imag();
    }
    matplot::plot(plotting_vector);
    matplot::hold(false);
    matplot::show();


    //signal1.show(PlottingKind::Spectre);
    //signal2.show(PlottingKind::Spectre);

    for (int i = 0; i < 500; i++) {
        auto freq = NP_DSP::ONE_D::UTILITY_MATH::getFreqByIdx(500, i);
        IC(freq, 1.0/freq, i);
    }
    //IC(*mask.base->vec);
    return 0;
}