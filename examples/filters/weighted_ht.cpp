#include <utility_math.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <complex>
#include <cstdlib>
#include <filters.hpp>
#include <utility_math.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;
    std::vector<std::complex<double>> resp;
    std::vector<std::complex<double>> filter;
    std::vector<double> plotting_vector;
    SignalT weights;

    for (int i = 0; i < 500; i++) {
        signal1.base->vec->push_back(std::rand());
        signal3.base->vec->push_back(0.0);
        plotting_vector.push_back(0.0);
        filter.push_back({0.0, 0.0});
        resp.push_back({0.0, 0.0});
    }

    double avg = 0.0;
    for (int i = 0; i < 500; i++){
        avg += signal1[i];
    }
    avg /= 500;
    for (int i = 0; i < 500; i++){
        signal1[i] -= avg;
    }
    
    for (int i = 0; i < 50; i++) {
        //signal2.base->vec->push_back(1.0);
        //resp.push_back({1.0, 0.0});
        weights.base->vec->push_back(0.0);
    }
    for (int i = 50; i < 450; i++){
        //resp.push_back({0.0, 0.0});
        weights.base->vec->push_back(1.0);
    }
    for (int i = 450; i < 500; i++) {
        //signal2.base->vec->push_back(1.0);
        //resp.push_back({1.0, 0.0});
        weights.base->vec->push_back(0.0);
    }

    NP_DSP::ONE_D::UTILITY_MATH::WeightedHilbertTransform
        <decltype(signal1), decltype(signal3), decltype(weights), NP_DSP::ONE_D::UTILITY_MATH::HTKind::Mull>
        (signal1, signal3, resp, filter, weights);

    matplot::plot(*signal1.base->vec);
    matplot::hold(true);
    for (int i = 0; i < 500; i++){
        plotting_vector[i] = filter[i].imag();
    }
    matplot::plot(*signal3.base->vec);
    matplot::hold(false);
    matplot::show();

    for (int i = 0; i < 500; i++) {
        auto freq = NP_DSP::ONE_D::UTILITY_MATH::getFreqByIdx(500, i);
        IC(freq, 1.0/freq, i);
    }

    return 0;
}