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

    NP_DSP::ONE_D::Tokenizers::instFreqNormSincReqTokenizer tokenizer;
    tokenizer.locality_coeff = 2.5;
    tokenizer.period_muller = 1.05;
    tokenizer.max_iter_number_for_filter = 3;
    tokenizer.debug = false;


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

