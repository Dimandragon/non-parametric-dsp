#pragma once

#include <cstddef>
#include <npdsp_concepts.hpp>
#include <signals.hpp>
#include <vector>
#include <inst_ampl_computers.hpp>
#include <phase_computers.hpp>
#include <inst_freq_computers.hpp>
#include <filters.hpp>
#include <integrators.hpp>
#include <derivators.hpp>

struct Token{
    size_t t;
    size_t mode_num;
    double val;
    double inst_freq;
    double inst_ampl;
    double phase;

    Token(){}
};

namespace NP_DSP::ONE_D::Tokenizers {
    struct InstFreqNormSincTokenizer
    {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        DataType data;
        DataType data_buffer;
        DataType compute_buffer;
        DataType compute_buffer2;
        DataType mode;
        DataType inst_freq;
        DataType inst_ampl;
        DataType phase;

        size_t modes_count;

        std::vector<double> freq_conv;
        std::vector<double> freq_conv_image;

        std::vector<Token> tokens;

        double period_muller = 1.0;
        double locality_coeff = 5.0;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
                <double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
                phase_computer_simple;

        INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_freq_computer =
                INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                        decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                        (integrator, derivator);

        INST_AMPL_COMPUTERS::HilbertTransformBased
                <UTILITY_MATH::HTKind::Mull> inst_ampl_computer;

        FILTERS::MonoFreqFilters<double, FILTERS::MonoInstFreqFilteringType::SincPaddedFIR> filter;


        template<typename DataT>
        void compute(const DataT & data_in){
            filter.gaussian_width_muller = locality_coeff;
            size_t iter_number = 0;

            if(data.size() != data_in.size()){
                data.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    data.base->vec->push_back(data_in[i]);
                }
            }
            if(data_buffer.size() != data_in.size()){
                data_buffer.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    data_buffer.base->vec->push_back(data_in[i]);
                }
            }
            if(compute_buffer.size() != data_in.size()){
                compute_buffer.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    compute_buffer.base->vec->push_back(0.0);
                }
            }
            if(compute_buffer2.size() != data_in.size()){
                compute_buffer2.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    compute_buffer2.base->vec->push_back(0.0);
                }
            }
            if(freq_conv.size() != data.size()){
                freq_conv.clear();
                for (int i = 0; i < data.size(); i++) {
                    freq_conv.push_back(1.0);
                }
            }
            if(freq_conv_image.size() != data.size()){
                freq_conv_image.clear();
                for (int i = 0; i < data.size(); i++) {
                    freq_conv_image.push_back(1.0);
                }
            }

            if(mode.size() != data.size()){
                mode.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    mode.base->vec->push_back(0.0);
                }
            }

            if(inst_freq.size() != data.size()){
                inst_freq.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    inst_freq.base->vec->push_back(0.0);
                }
            }

            if(inst_ampl.size() != data.size()){
                inst_ampl.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    inst_ampl.base->vec->push_back(0.0);
                }
            }

            if(phase.size() != data.size()){
                phase.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    phase.base->vec->push_back(0.0);
                }
            }

            tokens.clear();

            while(true){
                phase_computer_simple.compute(data, phase, nullptr);

                if(phase[data.size() - 1] > 6.28){
                    inst_freq_computer.compute(phase, inst_freq, nullptr);
                    double base_inst_freq = INST_FREQ_COMPUTERS::instFreqNorm(data, data_buffer, inst_freq, freq_conv, freq_conv_image);

                    phase_computer_simple.compute(data_buffer, phase, nullptr);

                    base_inst_freq = 1.0 /
                                     (static_cast<double>(data.size()) /
                                      (phase[data.size() - 1] / 2.0 / std::numbers::pi)) / period_muller;

                    for(int i = 0; i < data.size(); i++){
                        data[i] = data_buffer[i];
                    }

                    filter.freq = base_inst_freq;
                    filter.is_low_pass = true;
                    filter.compute(data, data_buffer, data_buffer);

                    for(int i = 0; i < data.size(); i++){
                        auto swap = data_buffer[i];
                        data_buffer[i] = data[i] - data_buffer[i];
                        data[i] = swap;
                    }

                    INST_FREQ_COMPUTERS::backInstFreqNorm(data_buffer, mode, freq_conv);

                    phase_computer_simple.compute(mode, phase, nullptr);

                    inst_freq_computer.compute(phase, inst_freq, nullptr);

                    inst_ampl_computer.compute(mode,  inst_ampl, nullptr);

                    auto i_temp = 0;
                    Token token_temp;
                    token_temp.mode_num = iter_number;
                    token_temp.inst_ampl = inst_ampl[i_temp];
                    token_temp.inst_freq = inst_freq[i_temp];
                    token_temp.phase = phase[i_temp];
                    token_temp.val = mode[i_temp];
                    token_temp.t = i_temp;
                    tokens.push_back(token_temp);

                    double current_mode = 0.0;
                    for (size_t i = 0; i < data.size(); i++){
                        if (phase[i] > current_mode + std::numbers::pi){
                            Token token;
                            token.mode_num = iter_number;
                            token.inst_ampl = inst_ampl[i];
                            token.inst_freq = inst_freq[i];
                            token.phase = phase[i];
                            token.val = mode[i];
                            token.t = i;
                            current_mode += std::numbers::pi;
                            tokens.push_back(token);
                        }
                    }

                    iter_number++;
                }
                else{
                    INST_FREQ_COMPUTERS::backInstFreqNorm(data, mode, freq_conv);

                    //std::cout << "last mode back inst freq norm " << iter_number << std::endl;
                    //modes[iter_number]->show(PlottingKind::Simple);
                    phase_computer_simple.compute(mode, phase, nullptr);
                    inst_freq_computer.compute(phase, inst_freq, nullptr);
                    inst_ampl_computer.compute(mode,  inst_ampl, nullptr);

                    auto i_temp = 0;
                    Token token_temp;
                    token_temp.mode_num = iter_number;
                    token_temp.inst_ampl = inst_ampl[i_temp];
                    token_temp.inst_freq = inst_freq[i_temp];
                    token_temp.phase = phase[i_temp];
                    token_temp.val = mode[i_temp];
                    token_temp.t = i_temp;
                    tokens.push_back(token_temp);

                    double current_phase = 0.0;
                    for (size_t i = 0; i < data.size(); i++){
                        if (phase[i] >= current_phase + std::numbers::pi){
                            Token token;
                            token.mode_num = iter_number;
                            token.inst_ampl = inst_ampl[i];
                            token.inst_freq = inst_freq[i];
                            token.phase = phase[i];
                            token.val = mode[i];
                            token.t = i;
                            current_phase += std::numbers::pi;
                            tokens.push_back(token);
                        }
                    }

                    iter_number++;
                    return;
                }
            }
        }

        std::vector<Token> getTokens(){
            return tokens;
        }
    };

    struct InstFreqNormSincReqTokenizer
    {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        DataType data;
        DataType data_buffer;
        DataType compute_buffer;
        DataType compute_buffer2;
        DataType mode;
        DataType inst_freq;
        DataType inst_ampl;
        DataType phase;
        std::vector<double> freq_conv;
        std::vector<double> freq_conv_image;

        double period_muller = 1.0;
        double locality_coeff = 5.0;
        double max_iter_number_for_filter = 10;

        bool debug = false;

        std::vector<Token> tokens;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
                <double, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
                phase_computer_der_atan;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
                <double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
                phase_computer_simple;

        INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage>
                inst_freq_computer =
                INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                        decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage>
                        (integrator, derivator);

        INST_AMPL_COMPUTERS::HilbertTransformBased
                <UTILITY_MATH::HTKind::Mull> inst_ampl_computer;

        FILTERS::RecursiveFilter<double, FILTERS::LocalFilteringType::SincResampled> filter;

        template<typename DataT>
        void compute(const DataT & data_in){
            filter.locality_coeff = locality_coeff;
            filter.period_muller = period_muller;
            //filter.inst_freq_computer = &inst_freq_computer;
            //filter.phase_computer = &phase_computer_simple;
            filter.debug = false;
            filter.max_iters = max_iter_number_for_filter;

            size_t iter_number = 0;

            DataType non_resampled_data;

            if(data.size() != data_in.size()){
                data.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    data.base->vec->push_back(data_in[i]);
                }
            }
            if(data_buffer.size() != data_in.size()){
                data_buffer.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    data_buffer.base->vec->push_back(data_in[i]);
                }
            }
            if(compute_buffer.size() != data_in.size()){
                compute_buffer.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    compute_buffer.base->vec->push_back(0.0);
                }
            }
            if(compute_buffer2.size() != data_in.size()){
                compute_buffer2.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    compute_buffer2.base->vec->push_back(0.0);
                }
            }
            if(freq_conv.size() != data.size()){
                freq_conv.clear();
                for (int i = 0; i < data.size(); i++) {
                    freq_conv.push_back(1.0);
                }
            }
            if(freq_conv_image.size() != data.size()){
                freq_conv_image.clear();
                for (int i = 0; i < data.size(); i++) {
                    freq_conv_image.push_back(1.0);
                }
            }

            if(mode.size() != data.size()){
                mode.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    mode.base->vec->push_back(0.0);
                }
            }

            if(inst_freq.size() != data.size()){
                inst_freq.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    inst_freq.base->vec->push_back(0.0);
                }
            }

            if(inst_ampl.size() != data.size()){
                inst_ampl.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    inst_ampl.base->vec->push_back(0.0);
                }
            }

            if(phase.size() != data.size()){
                phase.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    phase.base->vec->push_back(0.0);
                }
            }

            tokens.clear();

            while(true){
                phase_computer_simple.compute(data, phase, nullptr);

                if (debug){
                    std::cout << "compute first phase  " << iter_number << std::endl;
                    phase.show(PlottingKind::Simple);
                }

                if(phase[data.size() - 1] > 6.28){
                    filter.compute(data, data_buffer, &compute_buffer);

                    if (debug){
                        std::cout << "get low freq part  " << iter_number << std::endl;
                        data_buffer.show(PlottingKind::Simple);//, label.str());
                    }

                    for(int i = 0; i < data.size(); i++){
                        mode[i] = data[i] - data_buffer[i];
                        data[i] = data_buffer[i]; //data is filtered signal
                                        //data_buffer is mode
                    }
                    if (debug){
                        std::cout << "get mode " << iter_number << std::endl;
                        mode.show(PlottingKind::Simple);
                    }
                    if (debug){
                        std::cout << "get low_freq_part " << iter_number << std::endl;
                        data.show(PlottingKind::Simple);
                    }

                    phase_computer_simple.compute(mode, phase, nullptr);

                    inst_freq_computer.compute(phase, inst_freq, nullptr);

                    inst_ampl_computer.compute(mode, inst_ampl, nullptr);

                    auto i_temp = 0;
                    Token token_temp;
                    token_temp.mode_num = iter_number;
                    token_temp.inst_ampl = inst_ampl[i_temp];
                    token_temp.inst_freq = inst_freq[i_temp];
                    token_temp.phase = phase[i_temp];
                    token_temp.val = mode[i_temp];
                    token_temp.t = i_temp;
                    tokens.push_back(token_temp);

                    double current_mode = 0.0;
                    for (size_t i = 0; i < data.size(); i++){
                        if (phase[i] > current_mode + std::numbers::pi){
                            Token token;
                            token.mode_num = iter_number;
                            token.inst_ampl = inst_ampl[i];
                            token.inst_freq = inst_freq[i];
                            token.phase = phase[i];
                            token.val = mode[i];
                            token.t = i;
                            current_mode += std::numbers::pi;
                            tokens.push_back(token);
                        }
                    }

                    iter_number++;
                }
                else{
                    for (auto i = 0; i < data.size(); i++){
                        mode[i] = data[i];
                    }

                    phase_computer_simple.compute(mode, phase, nullptr);
                    inst_freq_computer.compute(phase, inst_freq, nullptr);
                    inst_ampl_computer.compute(mode, inst_ampl, nullptr);

                    auto i_temp = 0;
                    Token token_temp;
                    token_temp.mode_num = iter_number;
                    token_temp.inst_ampl = inst_ampl[i_temp];
                    token_temp.inst_freq = inst_freq[i_temp];
                    token_temp.phase = phase[i_temp];
                    token_temp.val = mode[i_temp];
                    token_temp.t = i_temp;
                    tokens.push_back(token_temp);

                    double current_phase = 0.0;
                    for (size_t i = 0; i < data.size(); i++){
                        if (phase[i] >= current_phase + std::numbers::pi){
                            Token token;
                            token.mode_num = iter_number;
                            token.inst_ampl = inst_ampl[i];
                            token.inst_freq = inst_freq[i];
                            token.phase = phase[i];
                            token.val = mode[i];
                            token.t = i;
                            current_phase += std::numbers::pi;
                            tokens.push_back(token);
                        }
                    }

                    iter_number++;
                    return;
                }
            }
        }

        std::vector<Token> getTokens(){
            return tokens;
        }
    };

    struct MakimaModeDecompositionBasedTokenizer
    {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        DataType data;
        DataType data_buffer;
        DataType compute_buffer;
        DataType compute_buffer2;
        DataType mode;
        DataType inst_freq;
        DataType inst_ampl;
        DataType phase;
        std::vector<double> freq_conv;
        std::vector<double> freq_conv_image;

        double max_iter_number_for_filter = 10;

        std::vector<double> phase_shifts 
            { 0.0 * std::numbers::pi, 0.1 * std::numbers::pi, 0.2 * std::numbers::pi, 
            0.3 * std::numbers::pi, 0.4 * std::numbers::pi, 0.5 * std::numbers::pi,
            0.6 * std::numbers::pi, 0.7 * std::numbers::pi, 0.8 * std::numbers::pi, 
            0.9 * std::numbers::pi};

        bool debug = false;

        std::vector<Token> tokens;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
                <double, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
                phase_computer_der_atan;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
                <double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
                phase_computer_simple;

        INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage>
                inst_freq_computer =
                INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                        decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage>
                        (integrator, derivator);

        INST_AMPL_COMPUTERS::HilbertTransformBased
                <UTILITY_MATH::HTKind::Mull> inst_ampl_computer;

        FILTERS::RecursiveFilter<double, FILTERS::LocalFilteringType::MakimaInterpolationExtremums> filter;

        template<typename DataT>
        void compute(const DataT & data_in){
            filter.filter.phase_shifts = phase_shifts;
            //filter.inst_freq_computer = &inst_freq_computer;
            //filter.phase_computer = &phase_computer_simple;
            filter.debug = false;
            filter.max_iters = max_iter_number_for_filter;

            size_t iter_number = 0;

            DataType non_resampled_data;

            if(data.size() != data_in.size()){
                data.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    data.base->vec->push_back(data_in[i]);
                }
            }
            if(data_buffer.size() != data_in.size()){
                data_buffer.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    data_buffer.base->vec->push_back(data_in[i]);
                }
            }
            if(compute_buffer.size() != data_in.size()){
                compute_buffer.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    compute_buffer.base->vec->push_back(0.0);
                }
            }
            if(compute_buffer2.size() != data_in.size()){
                compute_buffer2.base->vec->clear();
                for (int i = 0; i < data_in.size(); i++){
                    compute_buffer2.base->vec->push_back(0.0);
                }
            }
            if(freq_conv.size() != data.size()){
                freq_conv.clear();
                for (int i = 0; i < data.size(); i++) {
                    freq_conv.push_back(1.0);
                }
            }
            if(freq_conv_image.size() != data.size()){
                freq_conv_image.clear();
                for (int i = 0; i < data.size(); i++) {
                    freq_conv_image.push_back(1.0);
                }
            }

            if(mode.size() != data.size()){
                mode.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    mode.base->vec->push_back(0.0);
                }
            }

            if(inst_freq.size() != data.size()){
                inst_freq.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    inst_freq.base->vec->push_back(0.0);
                }
            }

            if(inst_ampl.size() != data.size()){
                inst_ampl.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    inst_ampl.base->vec->push_back(0.0);
                }
            }

            if(phase.size() != data.size()){
                phase.base->vec->clear();
                for (int i = 0; i < data.size(); i++) {
                    phase.base->vec->push_back(0.0);
                }
            }

            tokens.clear();

            while(true){
                phase_computer_simple.compute(data, phase, nullptr);


                if(phase[data.size() - 1] > 6.28){
                    filter.compute(data, data_buffer, &compute_buffer);

                    for(int i = 0; i < data.size(); i++){
                        mode[i] = data[i] - data_buffer[i];
                        data[i] = data_buffer[i]; //data is filtered signal
                                        //data_buffer is mode
                    }

                    phase_computer_simple.compute(mode, phase, nullptr);

                    inst_freq_computer.compute(phase, inst_freq, nullptr);

                    inst_ampl_computer.compute(mode, inst_ampl, nullptr);

                    auto i_temp = 0;
                    Token token_temp;
                    token_temp.mode_num = iter_number;
                    token_temp.inst_ampl = inst_ampl[i_temp];
                    token_temp.inst_freq = inst_freq[i_temp];
                    token_temp.phase = phase[i_temp];
                    token_temp.val = mode[i_temp];
                    token_temp.t = i_temp;
                    tokens.push_back(token_temp);

                    double current_mode = 0.0;
                    for (size_t i = 0; i < data.size(); i++){
                        if (phase[i] > current_mode + std::numbers::pi){
                            Token token;
                            token.mode_num = iter_number;
                            token.inst_ampl = inst_ampl[i];
                            token.inst_freq = inst_freq[i];
                            token.phase = phase[i];
                            token.val = mode[i];
                            token.t = i;
                            current_mode += std::numbers::pi;
                            tokens.push_back(token);
                        }
                    }

                    iter_number++;
                }
                else{
                    for (auto i = 0; i < data.size(); i++){
                        mode[i] = data[i];
                    }

                    phase_computer_simple.compute(mode, phase, nullptr);
                    inst_freq_computer.compute(phase, inst_freq, nullptr);
                    inst_ampl_computer.compute(mode, inst_ampl, nullptr);

                    auto i_temp = 0;
                    Token token_temp;
                    token_temp.mode_num = iter_number;
                    token_temp.inst_ampl = inst_ampl[i_temp];
                    token_temp.inst_freq = inst_freq[i_temp];
                    token_temp.phase = phase[i_temp];
                    token_temp.val = mode[i_temp];
                    token_temp.t = i_temp;
                    tokens.push_back(token_temp);

                    double current_phase = 0.0;
                    for (size_t i = 0; i < data.size(); i++){
                        if (phase[i] >= current_phase + std::numbers::pi){
                            Token token;
                            token.mode_num = iter_number;
                            token.inst_ampl = inst_ampl[i];
                            token.inst_freq = inst_freq[i];
                            token.phase = phase[i];
                            token.val = mode[i];
                            token.t = i;
                            current_phase += std::numbers::pi;
                            tokens.push_back(token);
                        }
                    }

                    iter_number++;
                    return;
                }
            }
        }

        std::vector<Token> getTokens(){
            return tokens;
        }
    };
}