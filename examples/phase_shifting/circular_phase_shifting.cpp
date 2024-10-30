#include "npdsp_concepts.hpp"
#include <signals.hpp>
#include <phase_shifters.hpp>
#include <derivators.hpp>
#include <cmath>
#include <utility_math.hpp>

int main(){
    using namespace NP_DSP::ONE_D;
    int size = 50;
    double res_ratio = 20;

    PHASE_SHIFTERS::HTBased ht_based;
    DERIVATORS::FTBased<DERIVATORS::FTDerivativeKind::Naive> ft_based1;
    DERIVATORS::FTBased<DERIVATORS::FTDerivativeKind::Rieze> ft_based2;
    DERIVATORS::FTBased<DERIVATORS::FTDerivativeKind::Weyl> ft_based3;

    PHASE_SHIFTERS::FracDiffsBasedSimple<decltype(ft_based1)> phase_shifter1;
    phase_shifter1.derivator = &ft_based1;
    PHASE_SHIFTERS::WithOversampling<decltype(phase_shifter1)> phase_shifter2;
    phase_shifter2.phase_shifter = &phase_shifter1;
    phase_shifter2.oversampling_ratio = 15;


    

    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2, signal_ref;

    for (int i = 0; i < size; i++){
        signal_ref.base->vec->push_back(std::rand());
    }

    size *= res_ratio;

    for (int i = 0; i < size; i++){
        signal1.base->vec->push_back(signal_ref.interpolate(double(i) /res_ratio, NP_DSP::ONE_D::SignalKind::Universal));
        signal2.base->vec->push_back(0.0);
    }
    matplot::plot(*signal1.base->vec);
    matplot::hold(true);
    for(int i = 0; i < 20; i++){
        double phase_shift = std::numbers::pi / 10 * (i);
        ft_based1.power = phase_shift/std::numbers::pi * 2.0;
        IC(phase_shift, ft_based1.power);
        ft_based1.compute(signal1, signal2, nullptr);
        NP_DSP::ONE_D::UTILITY_MATH::normalizeSTD(signal1, signal2);
        
        matplot::plot(*signal2.base->vec);
    }
    matplot::hold(false);
    matplot::show();
}