#include <modes_extractors.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <cstdlib>
#include <mode_colleretion_tester.cpp>
#include <iostream>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> data;

    for(int i = 0; i < 500; i++) {
        data.base->vec->push_back(std::rand());
    }

    NP_DSP::ONE_D::MODES_EXTRACTORS::InstFreqNormSincExtractorReq extractor;
    extractor.locality_coeff = 5;
    extractor.period_muller = 1.05;
    extractor.max_iter_number_for_filter = 3;
    extractor.debug = false;
    //extractor.non_opt_filter.period_muller = 1.0;
    //NP_DSP::ONE_D::MODES_EXTRACTORS::ByIterStopFunc stop;
    //data.show(NP_DSP::ONE_D::PlottingKind::Simple);
    //extractor.load(data);
    
    //data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    extractor.compute(data);

    std::cout << extractor.getModesCount() << std::endl;
    return 0;
}

