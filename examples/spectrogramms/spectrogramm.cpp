#include "matplot/freestanding/axes_functions.h"
#include "matplot/freestanding/plot.h"
#include <icecream.hpp>
#include <modes_extractors.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <cstdlib>
#include <mode_colleretion_tester.hpp>
#include <vector>
#include <matplot/matplot.h>
#include <iostream>
#include <math.h>
#include <spectrogramm.hpp>

bool save = true;

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> data;

    auto res_coeff = 2.0;

    for(int i = 0; i < 1000; i++) {
        double idx = i;
        data.base->vec->push_back(50 + std::sin(idx / res_coeff) * 50 + std::sin(idx * 3.14 / res_coeff) * 50 + std::cos(idx / 1000.0 * 3.14 / res_coeff) * 20
            + std::cos(idx * idx / 100000.0 / res_coeff) + std::cos(idx * idx / 100000.0 / res_coeff) * 200);
    }

    IC(*data.base->vec);
    if (save){
        data.show(NP_DSP::ONE_D::PlottingKind::Simple, 
            "/home/dmitry/projects/non-parametric-dsp/examples/spectrogramms/signal.svg");
    }
    else{
        data.show(NP_DSP::ONE_D::PlottingKind::Simple);
    }
    

    NP_DSP::ONE_D::MODES_EXTRACTORS::instFreqNormSincExtractorReq extractor;
    extractor.locality_coeff = 3;
    extractor.period_muller = 1.5;
    extractor.max_iter_number_for_filter = 5;
    extractor.debug = true;

    extractor.compute(data);

    for(int i = 0; i < extractor.modes.size(); i++) {
        std::stringstream path;
        path << "/home/dmitry/projects/non-parametric-dsp/examples/spectrogramms/mode" << i << ".png";

        if (save){
            extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple, path.str());
        }
        else{
            extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        }
    }

    matplot::plot(*extractor.inst_freqs[0]->base->vec);
    matplot::hold(true);
    
    for (int i = 1; i < extractor.getModesCount(); i++){
        matplot::plot(*extractor.inst_freqs[i]->base->vec);
    }
    matplot::hold(false);
    matplot::show();

    IC(extractor.modes.size());

    spectrogramm(extractor, 1000, 200);

    return 0;
}