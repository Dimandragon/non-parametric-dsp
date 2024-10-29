#include "npdsp_concepts.hpp"
#include <signals.hpp>
#include <phase_shifters.hpp>
#include <derivators.hpp>
#include <cmath>
#include <utility_math.hpp>

int main(){
    using namespace NP_DSP::ONE_D;
    auto size = 1000;

    GenericSignal<SimpleVecWrapper<double>, true> data1, data2, out1, out2, out3, out4;
    PHASE_SHIFTERS::HTBased ht_based;
    DERIVATORS::FTBased<DERIVATORS::FTDerivativeKind::Naive> ft_based1;
    DERIVATORS::FTBased<DERIVATORS::FTDerivativeKind::Rieze> ft_based2;
    DERIVATORS::FTBased<DERIVATORS::FTDerivativeKind::Weyl> ft_based3;

    PHASE_SHIFTERS::FracDiffsBasedSimple<decltype(ft_based1)> phase_shifter1;
    phase_shifter1.derivator = &ft_based1;
    PHASE_SHIFTERS::WithOversampling<decltype(phase_shifter1)> phase_shifter2;
    phase_shifter2.phase_shifter = &phase_shifter1;
    phase_shifter2.oversampling_ratio = 15;


    for (int i = 0; i < size; i++){
        data1.base->vec->push_back(std::rand());
        data2.base->vec->push_back(std::sin(i / 20.) + std::cos (i / 100.) * 2);
        out1.base->vec->push_back(0.0);
        out2.base->vec->push_back(0.0);
        out3.base->vec->push_back(0.0);
        out4.base->vec->push_back(0.0);
    }

    for(int i = 0; i < 200; i++){
        double phase_shift = std::numbers::pi / 100 * (i);
        //IC(phase_shift);
        ht_based.phase_shift = phase_shift;
        ft_based1.power = phase_shift/std::numbers::pi * 2.0;
        ft_based2.power = phase_shift/std::numbers::pi * 2.0;
        ft_based3.power = phase_shift/std::numbers::pi * 2.0;
        phase_shifter1.phase_shift = phase_shift;
        phase_shifter2.phase_shift = phase_shift;
        //ht_based.compute(data2, out1);
        ft_based1.compute(data1, out1, nullptr);
        phase_shifter1.compute(data1, out2, nullptr);
        phase_shifter2.compute(data1, out3, nullptr);
        ft_based2.compute(data1, out4, nullptr);
        //ft_based2.compute(data1, out2, nullptr);
        //ft_based3.compute(data1, out3, nullptr);
        NP_DSP::ONE_D::UTILITY_MATH::normalizeSTD(data1, out1);
        NP_DSP::ONE_D::UTILITY_MATH::normalizeSTD(data1, out2);
        NP_DSP::ONE_D::UTILITY_MATH::normalizeSTD(data1, out3);
        NP_DSP::ONE_D::UTILITY_MATH::normalizeSTD(data1, out4);
        /*data1.show(PlottingKind::Simple);
        out1.show(PlottingKind::Simple);
        out2.show(PlottingKind::Simple);*/
        matplot::plot(*data1.base->vec);
        matplot::hold(true);
        matplot::plot(*out1.base->vec);
        matplot::plot(*out2.base->vec);
        matplot::plot(*out3.base->vec);
        matplot::plot(*out4.base->vec);
        matplot::hold(false);
        matplot::show();
    }
}