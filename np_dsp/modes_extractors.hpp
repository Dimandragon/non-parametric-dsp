#pragma once
#include <npdsp_concepts.hpp>
#include <signals.hpp>
#include <vector>
#include <inst_ampl_computers.hpp>
#include <phase_computers.hpp>
#include <inst_freq_computers.hpp>
#include <filters.hpp>
#include <integrators.hpp>
#include <derivators.hpp>

namespace NP_DSP::ONE_D::MODES_EXTRACTORS {
    struct InstFreqNormSincExtractor
    {
        enum class ResamplingType{
            BackOnlyOnOut,
            FullBackAfterIter,
            BackForModeAfterIter,
            BackForSignalAfterIter,
        };

        ResamplingType resampling_type = ResamplingType::BackOnlyOnOut;

        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        DataType data;
        DataType data_buffer;
        DataType compute_buffer;
        DataType compute_buffer2;
        std::vector<DataType *> modes;
        std::vector<DataType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;
        std::vector<double> freq_conv;
        std::vector<double> freq_conv_image;

        double period_muller = 1.0;
        double locality_coeff = 5.0;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
                <double, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
                phase_computer_der_atan;

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

        FILTERS::NonOptPeriodBasedFilter<double,
                FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>
                non_opt_filter = FILTERS::NonOptPeriodBasedFilter<double,
                FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>(integrator);

        FILTERS::SincResLocalFilter<double> filter;

        template<typename DataT>
        void compute(const DataT & data_in){
            filter.locality_coeff = locality_coeff;
            size_t iter_number = 0;

            DataType non_resampled_data;

            auto prepare_memory_ext = [&](){
                for (int i = 0; i < data_in.size(); i++){
                    non_resampled_data.base->vec->push_back(data_in[i]);
                }
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
            };

            prepare_memory_ext();

            auto prepare_memory_int = [&](){
                if (iter_number + 1 > modes.size()) {
                    auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    modes.push_back(modes_new);
                    for (int i = 0; i < data.size(); i++) {
                        modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number + 1 > inst_freqs.size()) {
                    auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    inst_freqs.push_back(inst_freqs_new);
                    for (int i = 0; i < data.size(); i++) {
                        inst_freqs[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number + 1 > inst_ampls.size()) {
                    auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    inst_ampls.push_back(inst_ampls_new);
                    for (int i = 0; i < data.size(); i++) {
                        inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number + 1 > phases.size()) {
                    auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    phases.push_back(phase_new);
                    for (int i = 0; i < data.size(); i++) {
                        phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (modes[iter_number]->size() != data.size()) {
                    modes[iter_number]->base->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (inst_freqs[iter_number]->size() != data.size()) {
                    inst_freqs[iter_number]->base->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        inst_freqs[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (inst_ampls[iter_number]->size() != data.size()) {
                    inst_ampls[iter_number]->base->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (phases[iter_number]->size() != data.size()) {
                    phases[iter_number]->base->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
            };

            while(true){
                prepare_memory_int();
                phase_computer_simple.compute(data, *phases[iter_number], nullptr);

                //std::cout << "compute first phase  " << iter_number << std::endl;
                //phases[iter_number]->show(PlottingKind::Simple);

                if((*phases[iter_number])[data.size() - 1] > 6.28){
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);

                    double base_inst_freq = INST_FREQ_COMPUTERS::InstFreqNorm(data, data_buffer, *inst_freqs[iter_number], freq_conv, freq_conv_image);

                    phase_computer_simple.compute(data_buffer, *phases[iter_number], nullptr);

                    //std::cout << "compute resampled phase  " << iter_number << std::endl;
                    //phases[iter_number]->show(PlottingKind::Simple);

                    base_inst_freq = 1.0 /
                                     (static_cast<double>(data.size()) /
                                      ((*phases[iter_number])[data.size() - 1] / 2.0 / std::numbers::pi)) / period_muller;

                    //std::cout << "inst_freq_norm " << iter_number << std::endl;
                    //data_buffer.show(PlottingKind::Simple);

                    for(int i = 0; i < data.size(); i++){
                        data[i] = data_buffer[i];
                    }

                    filter.freq = base_inst_freq;
                    filter.is_low_pass = true;
                    filter.compute(data, data_buffer, nullptr);

                    //std::cout << "get low freq part  " << iter_number << std::endl;
                    //data_buffer.show(PlottingKind::Simple);//, label.str());
                    if (resampling_type == ResamplingType::BackForSignalAfterIter){
                        for (int i = 0; i < data.size(); i++){
                            (*modes[iter_number])[i] = data_buffer[i];
                        }
                        INST_FREQ_COMPUTERS::backInstFreqNorm(*modes[iter_number], data_buffer, freq_conv);
                        for (int i = 0; i < data.size(); i++){
                            (*modes[iter_number])[i] = data[i];
                        }
                        //INST_FREQ_COMPUTERS::backInstFreqNorm(*modes[iter_number], data, freq_conv);
                        for (int i = 0; i < data.size(); i++){
                            freq_conv[i] = 1.0;
                        }
                        for (int i = 0; i < data.size(); i++){
                            data[i] = non_resampled_data[i];
                        }
                    }

                    for(int i = 0; i < data.size(); i++){
                        //data_buffer is filtered signal, data is unfiltered
                        auto swap = data_buffer[i];
                        data_buffer[i] = data[i] - data_buffer[i];
                        data[i] = swap; //data is filtered signal
                                        //data_buffer is mode
                    }

                    //std::cout << "get mode " << iter_number << std::endl;
                    //data_buffer.show(PlottingKind::Simple);

                    if (resampling_type == ResamplingType::BackOnlyOnOut){
                        INST_FREQ_COMPUTERS::backInstFreqNorm(data_buffer, *modes[iter_number], freq_conv);
                    }
                    if (resampling_type == ResamplingType::FullBackAfterIter){
                        INST_FREQ_COMPUTERS::backInstFreqNorm(data_buffer, *modes[iter_number], freq_conv);
                        for (int i = 0; i < data.size(); i++){
                            data_buffer[i] = data[i];
                        }
                        INST_FREQ_COMPUTERS::backInstFreqNorm(data_buffer, data, freq_conv);
                        for (int i = 0; i < data.size(); i++){
                            freq_conv[i] = 1.0;
                        }
                    }
                    if (resampling_type == ResamplingType::BackForSignalAfterIter){
                        for (int i = 0; i < data.size(); i++){
                            (*modes[iter_number])[i] = data_buffer[i];
                        }
                        for (int i = 0; i < data.size(); i++){
                            non_resampled_data[i] = data[i];
                        }
                    }
                    if (resampling_type == ResamplingType::BackForModeAfterIter){
                        INST_FREQ_COMPUTERS::backInstFreqNorm(data_buffer, *modes[iter_number], freq_conv);
                        for (int i = 0; i < data.size(); i++){
                            freq_conv[i] = 1.0;
                        }
                        for(int i = 0; i < data.size(); i++){
                            non_resampled_data[i] -= (*modes[iter_number])[i];
                            data[i] = non_resampled_data[i];
                        }
                    }

                    //std::cout << "bask inst freq norm of mode " << iter_number << std::endl;
                    //modes[iter_number]->show(PlottingKind::Simple);

                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);

                    //std::cout << "compute result phase " << iter_number << std::endl;
                    //phases[iter_number]->show(PlottingKind::Simple);

                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);

                    //std::cout << "compute result inst_freq " << iter_number << std::endl;
                    //inst_freqs[iter_number]->show(PlottingKind::Simple);

                    inst_ampl_computer.compute(*modes[iter_number],  *inst_ampls[iter_number], nullptr);

                    //std::cout << "compute result inst_ampls " << iter_number << std::endl;
                    //inst_ampls[iter_number]->show(PlottingKind::Simple);

                    iter_number++;
                }
                else{
                    if (resampling_type == ResamplingType::BackOnlyOnOut){
                        INST_FREQ_COMPUTERS::backInstFreqNorm(data, *modes[iter_number], freq_conv);
                    }
                    else{
                        for (auto i = 0; i < data.size(); i++){
                            (*modes[iter_number])[i] = data[i];
                        }
                    }

                    //std::cout << "last mode back inst freq norm " << iter_number << std::endl;
                    //modes[iter_number]->show(PlottingKind::Simple);
                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    inst_ampl_computer.compute(*modes[iter_number],  *inst_ampls[iter_number], nullptr);

                    iter_number++;
                    return;
                }
            }
        }

        int getModesCount(){
            return static_cast<int>(modes.size());
        }
        int getDataSize(){
            return static_cast<int>(modes[0]->size());
        }
        std::vector<double> getMode(int idx){
            return *(modes[idx]->base->vec);
        }
        std::vector<double> getInstFreq(int idx){
            return *(inst_freqs[idx]->base->vec);
        }
        std::vector<double> getInstAmpl(int idx){
            return *(inst_ampls[idx]->base->vec);
        }
        std::vector<double> getPhase(int idx){
            return *(phases[idx]->base->vec);
        }
    };

    struct InstFreqNormSincExtractorReq
    {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        DataType data;
        DataType data_buffer;
        DataType compute_buffer;
        DataType compute_buffer2;
        std::vector<DataType *> modes;
        std::vector<DataType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;
        std::vector<double> freq_conv;
        std::vector<double> freq_conv_image;

        double period_muller = 1.0;
        double locality_coeff = 5.0;
        double max_iter_number_for_filter = 10;

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

        FILTERS::SincResLocalFilterWithResReq<double, decltype(phase_computer_simple), decltype(inst_freq_computer)> filter;

        template<typename DataT>
        void compute(const DataT & data_in){
            filter.locality_coeff = locality_coeff;
            filter.period_muller = period_muller;
            filter.inst_freq_computer = &inst_freq_computer;
            filter.phase_computer = &phase_computer_simple;
            filter.max_iters = max_iter_number_for_filter;

            size_t iter_number = 0;

            DataType non_resampled_data;

            auto prepare_memory_ext = [&](){
                for (int i = 0; i < data_in.size(); i++){
                    non_resampled_data.base->vec->push_back(data_in[i]);
                }
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
            };

            prepare_memory_ext();

            auto prepare_memory_int = [&](){
                if (iter_number + 1 > modes.size()) {
                    auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    modes.push_back(modes_new);
                    for (int i = 0; i < data.size(); i++) {
                        modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number + 1 > inst_freqs.size()) {
                    auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    inst_freqs.push_back(inst_freqs_new);
                    for (int i = 0; i < data.size(); i++) {
                        inst_freqs[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number + 1 > inst_ampls.size()) {
                    auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    inst_ampls.push_back(inst_ampls_new);
                    for (int i = 0; i < data.size(); i++) {
                        inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number + 1 > phases.size()) {
                    auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    phases.push_back(phase_new);
                    for (int i = 0; i < data.size(); i++) {
                        phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (modes[iter_number]->size() != data.size()) {
                    modes[iter_number]->base->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (inst_freqs[iter_number]->size() != data.size()) {
                    inst_freqs[iter_number]->base->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        inst_freqs[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (inst_ampls[iter_number]->size() != data.size()) {
                    inst_ampls[iter_number]->base->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (phases[iter_number]->size() != data.size()) {
                    phases[iter_number]->base->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
            };

            while(true){
                prepare_memory_int();
                phase_computer_simple.compute(data, *phases[iter_number], nullptr);

                if((*phases[iter_number])[data.size() - 1] > 6.28){
                    filter.compute(data, data_buffer, &compute_buffer);
                    for(int i = 0; i < data.size(); i++){
                        (*modes[iter_number])[i] = data[i] - data_buffer[i];
                        data[i] = data_buffer[i]; //data is filtered signal
                                        //data_buffer is mode
                    }

                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);

                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);

                    inst_ampl_computer.compute(*modes[iter_number],  *inst_ampls[iter_number], nullptr);

                    iter_number++;
                }
                else{
                    for (auto i = 0; i < data.size(); i++){
                        (*modes[iter_number])[i] = data[i];
                    }

                    //std::cout << "last mode back inst freq norm " << iter_number << std::endl;
                    //modes[iter_number]->show(PlottingKind::Simple);
                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    inst_ampl_computer.compute(*modes[iter_number],  *inst_ampls[iter_number], nullptr);

                    iter_number++;
                    return;
                }
            }
        }

        int getModesCount(){
            return static_cast<int>(modes.size());
        }
        int getDataSize(){
            return static_cast<int>(modes[0]->size());
        }
        std::vector<double> getMode(int idx){
            return *(modes[idx]->base->vec);
        }
        std::vector<double> getInstFreq(int idx){
            return *(inst_freqs[idx]->base->vec);
        }
        std::vector<double> getInstAmpl(int idx){
            return *(inst_ampls[idx]->base->vec);
        }
        std::vector<double> getPhase(int idx){
            return *(phases[idx]->base->vec);
        }
    };
}
