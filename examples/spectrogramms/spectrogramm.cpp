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

    for(int i = 0; i < 500; i++) {
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
    

    NP_DSP::ONE_D::MODES_EXTRACTORS::SOTAEMD extractor;
    //extractor.locality_coeff = 5;
    //extractor.period_muller = 1.2;
    extractor.max_iter_number_for_filter = 1;
    extractor.debug = false;
    extractor.oversampling_ratio_for_ft_der = 10.0;
    extractor.extremums_rotation_kind_e = NP_DSP::ONE_D::PHASE_SHIFTERS::RotateKind::Naive;
    extractor.phase_shifts = {0};
    for (int i = 0; i < 10; i++){
        extractor.phase_shifts.push_back(0.1 * i * std::numbers::pi);
    }
    extractor.interpolation_kind = NP_DSP::ONE_D::FILTERS::InterpolationKind::RBFMultiquadricAuto;

    //extractor.debug = true;
    extractor.compute(data);

    for(int i = 0; i < extractor.modes.size(); i++) {
        std::stringstream path;
        path << "/home/dmitry/projects/non-parametric-dsp/examples/spectrogramms/mode" << i << ".png";

        //if (save){
        //    extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple, path.str());
        //}
        //else{
            extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
            extractor.inst_freqs[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
            extractor.inst_ampls[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        // /}
    }

    /*matplot::plot(*extractor.inst_freqs[0]->base->vec);
    matplot::hold(true);
    
    for (int i = 1; i < extractor.getModesCount(); i++){
        matplot::plot(*extractor.inst_freqs[i]->base->vec);
    }
    matplot::hold(false);
    matplot::show();*/

    IC(extractor.modes.size());

    //spectrogramm(extractor, 500, 200);
    

    Spectrogramm spectrogramm;
    spectrogramm.setAxis(500, 500);
    spectrogramm.modifier = AmplModifier::log2;

    
    //spectrogramm.setBounds();

    spectrogramm.loadFromExtractor(extractor);
    spectrogramm.computeMatrix();
    // auto matrix1 = spectrogramm.matrix;
    spectrogramm.plot();
    //spectrogramm.smoothMatrix(20, 1);
    spectrogramm.smoothMatrixFast(20, 10);
    spectrogramm.plot();
    //IC(spectrogramm.computeRMSDistance(*matrix1), spectrogramm.computeL2Distance(*matrix1));

    return 0;
}