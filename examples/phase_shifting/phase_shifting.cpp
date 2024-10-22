#include "npdsp_concepts.hpp"
#include <signals.hpp>
#include <phase_shifters.hpp>
#include <derivators.hpp>
#include <cmath>

int main(){
    using namespace NP_DSP::ONE_D;
    auto size = 1000;

    GenericSignal<SimpleVecWrapper<double>, true> data1, data2, out1, out2;
    PHASE_SHIFTERS::HTBased ht_based;
    DERIVATORS::FTBased ft_based;

    for (int i = 0; i < size; i++){
        data1.base->vec->push_back(std::rand());
        data2.base->vec->push_back(std::sin(i / 20.) + std::cos (i / 100.) * 2);
        out1.base->vec->push_back(0.0);
        out2.base->vec->push_back(0.0);
    }

    for(int i = 0; i < 100; i++){
        double phase_shift = std::numbers::pi / 100 * (i) / 2.0;
        IC(phase_shift);
        ht_based.phase_shift = phase_shift;
        ft_based.power = phase_shift/std::numbers::pi;
        ht_based.compute(data2, out1);
        ft_based.compute(data2, out2, nullptr);
        /*data1.show(PlottingKind::Simple);
        out1.show(PlottingKind::Simple);
        out2.show(PlottingKind::Simple);*/
        matplot::plot(*data2.base->vec);
        matplot::hold(true);
        matplot::plot(*out1.base->vec);
        matplot::plot(*out2.base->vec);
        matplot::hold(false);
        matplot::show();
    }
}