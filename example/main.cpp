#include <modes_extractors.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <cstdlib>
#include <mode_colleretion_tester.hpp>
#include <iostream>
#include <tokenizer.hpp>

void tokenizerTest(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> data;
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> aaa;

    for(int i = 0; i < 500; i++) {
        data.base->vec->push_back(std::rand());
    }

    NP_DSP::ONE_D::Tokenizers::MakimaModeDecompositionBasedTokenizer tokenizer;
    tokenizer.phase_shifts = {};
    for (int i = 0; i < 100; i++){
        tokenizer.phase_shifts.push_back(0.01 * i * std::numbers::pi / 2.0);
    }
    tokenizer.max_iter_number_for_filter = 3;
    tokenizer.debug = false;


    tokenizer.compute(data);

    auto tokens = tokenizer.getTokens();

    for(const auto & token: tokens){
        std::cout << token.mode_num << " " << token.t << " " <<
            token.val << " " << token.inst_freq << " " <<
            token.inst_ampl << " " << token.phase << std::endl;  
    }
}

void extractorTest(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> data;

    for(int i = 0; i < 500; i++) {
        data.base->vec->push_back(std::rand());
    }

    NP_DSP::ONE_D::MODES_EXTRACTORS::MakimaBasedModeDecomposition extractor;
    extractor.max_iter_number_for_filter = 3;
    extractor.phase_shifts = {};
    for (int i = 0; i < 100; i++){
        extractor.phase_shifts.push_back(0.01 * i * std::numbers::pi);
    }
    extractor.debug = false;
    //extractor.non_opt_filter.period_muller = 1.0;
    //NP_DSP::ONE_D::MODES_EXTRACTORS::ByIterStopFunc stop;
    //data.show(NP_DSP::ONE_D::PlottingKind::Simple);
    //extractor.load(data);
    
    //data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    extractor.compute(data);

    std::cout << extractor.getModesCount() << std::endl;
}

int main(){
    tokenizerTest();
    extractorTest();
    return 0;
}

