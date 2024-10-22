#include <icecream.hpp>
#include <modes_extractors.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <cstdlib>
#include <mode_colleretion_tester.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> data;

    for(int i = 0; i < 500; i++) {
        data.base->vec->push_back(std::rand());
    }

    IC(*data.base->vec);
    data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::MODES_EXTRACTORS::instFreqNormSincExtractorReq extractor;
    extractor.locality_coeff = 5;
    extractor.period_muller = 1.15;
    extractor.max_iter_number_for_filter = 3;
    extractor.debug = false;
    //extractor.non_opt_filter.period_muller = 1.0;
    //NP_DSP::ONE_D::MODES_EXTRACTORS::ByIterStopFunc stop;
    //data.show(NP_DSP::ONE_D::PlottingKind::Simple);
    //extractor.load(data);
    
    //data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    extractor.compute(data);

    //NP_DSP::ONE_D::MODES_EXTRACTORS::computeReqDouble(extractor, stop);
    IC(extractor.modes.size());
    for(int i = 0; i < extractor.modes.size(); i++) {
        extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        //extractor.phases[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        //extractor.inst_freqs[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        //extractor.inst_ampls[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
    }
    return 0;
}

