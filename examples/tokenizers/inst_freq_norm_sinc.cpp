#include <modes_extractors.hpp>
//#include <cstdlib>
#include <iostream>
#include <tokenizer.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> data;
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> aaa;


    for(int i = 0; i < 500; i++) {
        data.base->vec->push_back(std::rand());
    }
    for(int i = 0; i < 5; i++) {
        aaa.base->vec->push_back(std::rand());
    }

    aaa.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::Tokenizers::InstFreqNormSincTokenizer tokenizer;
    tokenizer.locality_coeff = 2.5;
    tokenizer.period_muller = 1.0;

    tokenizer.compute(data);

    auto tokens = tokenizer.getTokens();

    //IC(tokens[0]);
    //IC(tokens);

    for(const auto & token: tokens){
        std::cout << token.mode_num << " " << token.t << " " <<
            token.val << " " << token.inst_freq << " " <<
            token.inst_ampl << " " << token.phase << std::endl;  
    }

    return 0;
}

