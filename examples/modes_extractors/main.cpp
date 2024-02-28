#include <icecream.hpp>

import modes_extractors;
import signals;
import npdsp_concepts;
import <vector>;
import <cstdlib>;

int main(){
    auto data = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});

    for(int i = 0; i < 50; i++) {
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(data.base)->vec->push_back(std::rand() / 10000);
        //data.base->vec->push_back(std::rand() / 100000.);
    }

    IC(*static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(data.base)->vec);

    NP_DSP::ONE_D::MODES_EXTRACTORS::MainExtractorDouble extractor;
    data.show(NP_DSP::ONE_D::PlottingKind::Simple);
    //extractor.load(data);//todo

    data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    extractor.compute();

    for(int i = 0; i < extractor.modes.size(); i++) {
        extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        extractor.phases[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        //extractor.inst_freqs[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        extractor.inst_ampls[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
    }

    return 0;
}