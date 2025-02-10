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

    NP_DSP::ONE_D::Tokenizers::SOTAEMDBasedTokenizer tokenizer;
    tokenizer.max_iter_number_for_filter = 5;
    tokenizer.debug = false;
    tokenizer.oversampling_ratio_for_ft_der = 10.0;
    tokenizer.extremums_rotation_kind_e = NP_DSP::ONE_D::PHASE_SHIFTERS::RotateKind::Naive;
    tokenizer.phase_shifts = {0};
    for (int i = 0; i < 100; i++){
        tokenizer.phase_shifts.push_back(0.01 * i * std::numbers::pi);
    }
    tokenizer.interpolation_kind_e = NP_DSP::ONE_D::FILTERS::InterpolationKind::RBFTPS;



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

    NP_DSP::ONE_D::MODES_EXTRACTORS::SOTAEMD extractor;
    extractor.max_iter_number_for_filter = 5;
    extractor.debug = false;
    extractor.oversampling_ratio_for_ft_der = 10.0;
    extractor.extremums_rotation_kind_e = NP_DSP::ONE_D::PHASE_SHIFTERS::RotateKind::Naive;
    extractor.phase_shifts = {0};
    for (int i = 0; i < 100; i++){
        extractor.phase_shifts.push_back(0.01 * i * std::numbers::pi);
    }
    extractor.interpolation_kind_e = NP_DSP::ONE_D::FILTERS::InterpolationKind::RBFTPS;

    extractor.compute(data);

    std::cout << extractor.getModesCount() << std::endl;
}

int main(){
    tokenizerTest();
    extractorTest();
    return 0;
}

