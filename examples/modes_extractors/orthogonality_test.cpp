#include <icecream.hpp>
#include <modes_extractors.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <cstdlib>
#include <mode_colleretion_tester.hpp>
#include <vector>
#include <matplot/matplot.h>
#include <cmath>
#include <iostream>
#include <math.h>

bool save = true;

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> data;

    for(int i = 0; i < 200; i++) {
        data.base->vec->push_back(std::rand() + std::sin(i / 500) * 2500000000. + std::sin(i / 20) * 100000000.);
    }

    IC(*data.base->vec);
    if (save){
        data.show(NP_DSP::ONE_D::PlottingKind::Simple, 
            "/home/dmitry/projects/non-parametric-dsp/examples/modes_extractors/signal.png");
    }
    else{
        data.show(NP_DSP::ONE_D::PlottingKind::Simple);
    }
    

    NP_DSP::ONE_D::MODES_EXTRACTORS::MakimaBasedModeDecomposition extractor;
    //extractor.locality_coeff = 5;
    //extractor.period_muller = 1.2;
    extractor.max_iter_number_for_filter = 5;
    extractor.debug = true;
    extractor.phase_shifts = {};
    for (int i = 0; i < 100; i++){
        extractor.phase_shifts.push_back(0.01 * i * std::numbers::pi);
    }

    extractor.compute(data);

    IC(extractor.modes.size());
    for(int i = 0; i < extractor.modes.size(); i++) {
        std::stringstream path;
        path << "/home/dmitry/projects/non-parametric-dsp/examples/modes_extractors/mode" << i << ".png";

        if (save){
            extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple, path.str());
        }
        else{
            extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        }
    }
    
    std::vector<std::vector<double>> orthogonality_distr;


    std::cout << computeOrthogonality<decltype(data), decltype(data)>
            (*(extractor.modes[0]), *(extractor.modes[1])) << " " 
            << computeOrthogonality<decltype(data), decltype(data)>
            (*(extractor.modes[1]), *(extractor.modes[0])) << std::endl;

    for (int i = 0; i < extractor.modes.size(); i++){
        orthogonality_distr.push_back({});
        for (int j = 0; j < extractor.modes.size(); j++){
            //std::cout << i << " " << j << std::endl;

            std::cout << computeOrthogonality<decltype(data), decltype(data)>
            (*(extractor.modes[i]), *(extractor.modes[j])) << " ";
            
            orthogonality_distr[i].push_back(std::round(100. * computeOrthogonality<decltype(data), decltype(data)>
            (*(extractor.modes[i]), *(extractor.modes[j]))) / 100.);
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    matplot::heatmap(orthogonality_distr);
    matplot::show();
    if (save){
        matplot::save("/home/dmitry/projects/non-parametric-dsp/examples/modes_extractors/orthogonality_test.png");
    }

    return 0;
}