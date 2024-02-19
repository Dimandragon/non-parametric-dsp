module;

import npdsp_concepts;
import signals;
import <vector>;
import inst_ampl_computers;
import phase_computers;
import inst_freq_computers;
import filters;
import integrators;
import derivators;

export module modes_extractors;




namespace NP_DSP{
    namespace ONE_D{
        namespace MODES_EXTRACTORS{
            export
            template <InstFreqComputer InstFreqComputerT>
            struct LinearExtractor{
                //using InstFreqComputerT
            };

            export 
            struct MainExtractor{
                using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
                DataType data;
                std::vector<DataType *> modes;
                std::vector<DataType *> inst_freqs;
                std::vector<DataType *> inst_ampls;
                std::vector<DataType *> phases;

                INTEGRATORS::Riman<DataType, DataType, INTEGRATORS::PolygonType::ByPoint> integrator;
                DERIVATORS::FinniteDifference<DataType, DataType, DERIVATORS::FinniteDifferenceType::Central> derivator;

                PHASE_COMPUTERS::ExtremumsBasedNonOpt<DataType, DataType> phase_computer;

                INST_FREQ_COMPUTERS::ComputedOnPhase<DataType, DataType, decltype(integrator),
                    decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                        inst_freq_computer;

                INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq<DataType, DataType, DataType,
                    DataType, decltype(integrator), decltype(derivator),
                        INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                            inst_ampl_computer;

                FILTERS::OptPeriodBasedFilter<DataType, DataType, DataType, FILTERS::FilteringType::DerivativeBased,
                    decltype(integrator), decltype(derivator), FILTERS::InstFreqComputerKind::phase_based_time_average,
                        FILTERS::PhaseComputingKind::extremums_based_non_opt> filter;

                ~MainExtractor(){
                    for (int i = 0; i < modes.size(); i++){
                        delete modes[i];
                    }
                    for (int i = 0; i < inst_freqs.size(); i++){
                        delete inst_freqs[i];
                    }
                    for (int i = 0; i < phases.size(); i++){
                        delete phases[i];
                    }
                    for (int i = 0; i < inst_ampls.size(); i++){
                        delete inst_ampls[i];
                    }
                }

                template<Signal DataInT>
                void load(const DataInT & data_in){
                    data.base->vec->clear();
                    for (auto i = 0; i < data_in.size(); i++){
                        data.base->vec->push_back(data.size());
                    }
                    for (int i = 0; i < modes.size(); i++){
                        modes[i]->base->vec->clear();
                    }
                    for (int i = 0; i < inst_freqs.size(); i++){
                        inst_freqs[i]->base->vec->clear();
                    }
                    for (int i = 0; i < phases.size(); i++){
                        phases[i]->base->vec->clear();
                    }
                    for (int i = 0; i < inst_ampls.size(); i++){
                        inst_ampls[i]->base->vec->clear();
                    }
                }

                bool stop(size_t iter_number){
                    if ((*phases[iter_number])[data.size()-1] < std::numbers::pi * 2){
                        return true;
                    }
                    else{
                        return false;
                    }
                }

                void compute(){
                    size_t iter_number = 0;
                    for(;;){
                        if (iter_number > modes.size() - 1) {
                            auto * modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                            modes.push_back(modes_new);
                            for (int i = 0; i < data.size(); i++) {
                                modes[iter_number]->base->vec->push_back(0.0);
                            }
                        }
                        if (iter_number > inst_freqs.size() - 1) {
                            auto * inst_freqs_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                            inst_freqs.push_back(inst_freqs_new);
                            for (int i = 0; i < data.size(); i++) {
                                inst_freqs[iter_number]->base->vec->push_back(0.0);
                            }
                        }
                        if (iter_number > inst_ampls.size() - 1) {
                            auto * inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                            inst_ampls.push_back(inst_ampls_new);
                            for (int i = 0; i < data.size(); i++) {
                                inst_ampls[iter_number]->base->vec->push_back(0.0);
                            }
                        }
                        if (iter_number > phases.size() - 1) {
                            auto * phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
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
                        GENERAL::Nil nil;

                        phase_computer.compute(data, *phases[iter_number], nil);

                        inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nil);

                        filter.compute(data, *modes[iter_number], *inst_freqs[iter_number]);

                        for (int i = 0; i < data.size(); i++) {
                            (*modes[iter_number])[i] = data[i] - (*modes[iter_number])[i];
                        }
                        for (int i = 0; i < data.size(); i++) {
                            data[i] = data[i] - (*modes[iter_number])[i];
                        }

                        DataType computer_buffer;
                        for (int i = 0; i < data.size(); i++) {
                            computer_buffer.base->vec->push_back(0.0);
                        }
                        phase_computer.compute(*modes[iter_number], *phases[iter_number], nil);
                        inst_ampl_computer.inst_freq = inst_freqs[iter_number];
                        inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], computer_buffer);

                        if (stop(iter_number)){
                            iter_number++;

                            if (iter_number > modes.size() - 1) {
                                auto * modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                                modes.push_back(modes_new);
                                for (int i = 0; i < data.size(); i++) {
                                    modes[iter_number]->base->vec->push_back(0.0);
                                }
                            }
                            if (iter_number > inst_freqs.size() - 1) {
                                auto * inst_freqs_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                                inst_freqs.push_back(inst_freqs_new);
                                for (int i = 0; i < data.size(); i++) {
                                    inst_freqs[iter_number]->base->vec->push_back(0.0);
                                }
                            }
                            if (iter_number > inst_ampls.size() - 1) {
                                auto * inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                                inst_ampls.push_back(inst_ampls_new);
                                for (int i = 0; i < data.size(); i++) {
                                    inst_ampls[iter_number]->base->vec->push_back(0.0);
                                }
                            }
                            if (iter_number > phases.size() - 1) {
                                auto * phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
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

                            for(int i = 0; i < data.size(); i++) {
                                (*modes[iter_number])[i] = data[i];
                            }

                            phase_computer.compute(*modes[iter_number], *phases[iter_number], nil);
                            inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nil);
                            inst_ampl_computer.inst_freq = inst_freqs[iter_number];
                            inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], computer_buffer);

                            break;
                        }
                        iter_number++;
                    }
                }

                //todo test
            };
        }
    }
}