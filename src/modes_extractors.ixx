module;

#include <icecream.hpp>

export module modes_extractors;

import npdsp_concepts;
import signals;
import <vector>;
import inst_ampl_computers;
import phase_computers;
import inst_freq_computers;
import filters;
import integrators;
import derivators;

namespace NP_DSP::ONE_D::MODES_EXTRACTORS {

    export
    struct MainExtractor {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<DataType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Central> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
            <double, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
                phase_computer;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
            <double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
                phase_computer_for_mode;

        INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
            decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_freq_computer =
                    INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                        decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                            (integrator, derivator);

        INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq<double, decltype(integrator), decltype(derivator),
            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_ampl_computer;


        /*FILTERS::OptPeriodBasedFilter<DataType, DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::DerArctg,
            FILTERS::FilteringType::DerivativeBased, decltype(integrator), decltype(derivator), FILTERS::InstFreqComputerKind::phase_based_time_average,
                FILTERS::PhaseComputingKind::extremums_based_non_opt>
                    filter = FILTERS::OptPeriodBasedFilter<DataType, DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::DerArctg, FILTERS::FilteringType::DerivativeBased,
                        decltype(integrator), decltype(derivator), FILTERS::InstFreqComputerKind::phase_based_time_average,
                            FILTERS::PhaseComputingKind::extremums_based_non_opt>(integrator, derivator);*/

        FILTERS::NonOptPeriodBasedFilter<double, FILTERS::FilteringType::Median,
            decltype(integrator), FILTERS::InstFreqKind::Average>
                filter1;

        FILTERS::NonOptPeriodBasedFilter<double, FILTERS::FilteringType::ValueBased,
            decltype(integrator), FILTERS::InstFreqKind::Average>
                filter2;

        ~MainExtractor() {
            for (int i = 0; i < modes.size(); i++) {
                delete modes[i];
            }
            for (int i = 0; i < inst_freqs.size(); i++) {
                delete inst_freqs[i];
            }
            for (int i = 0; i < phases.size(); i++) {
                delete phases[i];
            }
            for (int i = 0; i < inst_ampls.size(); i++) {
                delete inst_ampls[i];
            }
        }

        template<Signal DataInT>
        void load(const DataInT& data_in) {
            data.base->vec->clear();
            for (auto i = 0; i < data_in.size(); i++) {
                data.base->vec->push_back(data_in[i]);
            }
            for (int i = 0; i < modes.size(); i++) {
                modes[i]->base->vec->clear();
            }
            for (int i = 0; i < inst_freqs.size(); i++) {
                inst_freqs[i]->base->vec->clear();
            }
            for (int i = 0; i < phases.size(); i++) {
                phases[i]->base->vec->clear();
            }
            for (int i = 0; i < inst_ampls.size(); i++) {
                inst_ampls[i]->base->vec->clear();
            }
            for (int i = 0; i < inst_freqs.size(); i++) {
                inst_freqs[i]->base->vec->clear();
            }
        }

        bool stop(size_t iter_number) {
            IC(iter_number, (*phases[iter_number])[data.size()-1]);
            if ((*phases[iter_number])[data.size() - 1] < std::numbers::pi * 2.0) {
                return true;
            } else {
                return false;
            }
        }

        void compute() {
            size_t iter_number = 0;
            for (;;) {
                //data.show(PlottingKind::Simple);
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

                phase_computer.compute(data, *phases[iter_number], nullptr);
                //phases[iter_number]->show(PlottingKind::Simple);

                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                filter2.compute(data, *modes[iter_number], inst_freqs[iter_number]);


                //modes[iter_number]->show(PlottingKind::Simple);
                for (int i = 0; i < data.size(); i++) {
                    (*modes[iter_number])[i] = data[i] - (*modes[iter_number])[i];
                }
                for (int i = 0; i < data.size(); i++) {
                    data[i] = data[i] - (*modes[iter_number])[i];
                }
                //modes[iter_number]->show(PlottingKind::Simple);

                DataType computer_buffer;
                for (int i = 0; i < data.size(); i++) {
                    computer_buffer.base->vec->push_back(0.0);
                }
                phase_computer_for_mode.compute(*modes[iter_number], *phases[iter_number], nullptr);
                //phases[iter_number]->show(PlottingKind::Simple);
                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                inst_ampl_computer.inst_freq = inst_freqs[iter_number];
                inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);

                //inst_ampls[iter_number]->show(PlottingKind::Simple);

                if (stop(iter_number)) {
                    iter_number++;

                    if (iter_number > modes.size() - 1) {
                        auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                        modes.push_back(modes_new);
                        for (int i = 0; i < data.size(); i++) {
                            modes[iter_number]->base->vec->push_back(0.0);
                        }
                    }
                    if (iter_number > inst_freqs.size() - 1) {
                        auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                        inst_freqs.push_back(inst_freqs_new);
                        for (int i = 0; i < data.size(); i++) {
                            inst_freqs[iter_number]->base->vec->push_back(0.0);
                        }
                    }
                    if (iter_number > inst_ampls.size() - 1) {
                        auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                        inst_ampls.push_back(inst_ampls_new);
                        for (int i = 0; i < data.size(); i++) {
                            inst_ampls[iter_number]->base->vec->push_back(0.0);
                        }
                    }
                    if (iter_number > phases.size() - 1) {
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

                    for (int i = 0; i < data.size(); i++) {
                        (*modes[iter_number])[i] = data[i];
                    }

                    phase_computer.compute(*modes[iter_number], *phases[iter_number], nullptr);
                    //phases[iter_number]->show(PlottingKind::Simple);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    //inst_freqs[iter_number]->show(PlottingKind::Simple);
                    inst_ampl_computer.inst_freq = inst_freqs[iter_number];
                    inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);
                    //inst_ampls[iter_number]->show(PlottingKind::Simple);
                    //IC(iter_number, *phases[iter_number][]);
                    break;
                }
                iter_number++;
            }
        }
    };

    export
    struct MainExtractorDouble {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        using InstFreqType = GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<InstFreqType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Central> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
            <double, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
                phase_computer;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
            <double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
                phase_computer_for_mode;

        INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
            decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                inst_freq_computer =
                    INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                        decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                            (integrator, derivator);

        INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq<double, decltype(integrator), decltype(derivator),
            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                inst_ampl_computer;

  
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, decltype(phase_computer)>
                inst_freq_computer_for_opt_filter = NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
                    NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, 
                        decltype(phase_computer)> (integrator, derivator, phase_computer);

        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, decltype(phase_computer_for_mode)>
                inst_freq_computer_for_mode_for_opt_filter = NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), 
                    decltype(derivator), NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, 
                        decltype(phase_computer_for_mode)> (integrator, derivator, phase_computer_for_mode);
        
        NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter = NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);

        NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter2 = NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);
        
        NP_DSP::ONE_D::FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
            decltype(non_opt_filter2)> filter2 = NP_DSP::ONE_D::FILTERS::CascadeFilter
                <double, decltype(non_opt_filter), 
                    decltype(non_opt_filter2)>(non_opt_filter, non_opt_filter2);

        /*NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), decltype(inst_freq_computer_for_opt_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)> 
                filter2 = NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), 
                    decltype(inst_freq_computer_for_opt_filter), decltype(phase_computer), 
                        decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)> 
                            (cascade_filter, inst_freq_computer_for_opt_filter, phase_computer, 
                                inst_freq_computer_for_mode_for_opt_filter, phase_computer_for_mode);*/

        
        //FILTERS::CascadeFilter

        ~MainExtractorDouble() {
            for (int i = 0; i < modes.size(); i++) {
                delete modes[i];
            }
            for (int i = 0; i < inst_freqs.size(); i++) {
                delete inst_freqs[i];
            }
            for (int i = 0; i < phases.size(); i++) {
                delete phases[i];
            }
            for (int i = 0; i < inst_ampls.size(); i++) {
                delete inst_ampls[i];
            }
        }

        template<Signal DataInT>
        void load(const DataInT& data_in) {
            data.base->vec->clear();
            for (auto i = 0; i < data_in.size(); i++) {
                data.base->vec->push_back(data_in[i]);
            }
            for (int i = 0; i < modes.size(); i++) {
                modes[i]->base->vec->clear();
            }
            for (int i = 0; i < inst_freqs.size(); i++) {
                inst_freqs[i]->base->vec->clear();
            }
            for (int i = 0; i < phases.size(); i++) {
                phases[i]->base->vec->clear();
            }
            for (int i = 0; i < inst_ampls.size(); i++) {
                inst_ampls[i]->base->vec->clear();
            }
        }

        bool stop(size_t iter_number) {
            IC(iter_number, (*phases[iter_number])[data.size()-1]);
            if ((*phases[iter_number])[data.size() - 1] < std::numbers::pi * 2.0) {
                return true;
            } else {
                return false;
            }
        }

        void compute() {
            inst_freq_computer.variability = 0.5;
            size_t iter_number = 0;
            for (;;) {
                //data.show(PlottingKind::Simple);
                if (iter_number + 1 > modes.size()) {
                    auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    modes.push_back(modes_new);
                    for (int i = 0; i < data.size(); i++) {
                        modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number + 1 > inst_freqs.size()) {
                    auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                    inst_freqs.push_back(inst_freqs_new);
                    for (int i = 0; i < data.size(); i++) {
                        inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
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
                        inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
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
                auto size_data = data.size();
                auto out_size = phases[iter_number]->size();
                phase_computer.compute(data, *phases[iter_number], nullptr);
                //
                phases[iter_number]->show(PlottingKind::Simple);

                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                filter2.compute(data, *modes[iter_number], inst_freqs[iter_number]);

                //
                modes[iter_number]->show(PlottingKind::Simple);
                for (int i = 0; i < data.size(); i++) {
                    (*modes[iter_number])[i] = data[i] - (*modes[iter_number])[i];
                }
                for (int i = 0; i < data.size(); i++) {
                    data[i] = data[i] - (*modes[iter_number])[i];
                }
                //
                modes[iter_number]->show(PlottingKind::Simple);

                DataType computer_buffer;
                for (int i = 0; i < data.size(); i++) {
                    computer_buffer.base->vec->push_back(0.0);
                }

                phase_computer_for_mode.compute(*modes[iter_number], *phases[iter_number], nullptr);
                //phases[iter_number]->show(PlottingKind::Simple);
                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                inst_ampl_computer.inst_freq_double = inst_freqs[iter_number];
                inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);

                //inst_ampls[iter_number]->show(PlottingKind::Simple);

                if (stop(iter_number)) {
                    iter_number++;

                    if (iter_number > modes.size() - 1) {
                        auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                        modes.push_back(modes_new);
                        for (int i = 0; i < data.size(); i++) {
                            modes[iter_number]->base->vec->push_back(0.0);
                        }
                    }
                    if (iter_number > inst_freqs.size() - 1) {
                        auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                        inst_freqs.push_back(inst_freqs_new);
                        for (int i = 0; i < data.size(); i++) {
                            inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                        }
                    }
                    if (iter_number > inst_ampls.size() - 1) {
                        auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                        inst_ampls.push_back(inst_ampls_new);
                        for (int i = 0; i < data.size(); i++) {
                            inst_ampls[iter_number]->base->vec->push_back(0.0);
                        }
                    }
                    if (iter_number > phases.size() - 1) {
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
                            inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
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

                    for (int i = 0; i < data.size(); i++) {
                        (*modes[iter_number])[i] = data[i];
                    }

                    phase_computer.compute(*modes[iter_number], *phases[iter_number], nullptr);
                    //phases[iter_number]->show(PlottingKind::Simple);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    //inst_freqs[iter_number]->show(PlottingKind::Simple);
                    inst_ampl_computer.inst_freq_double = inst_freqs[iter_number];
                    inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);
                    //inst_ampls[iter_number]->show(PlottingKind::Simple);
                    //IC(iter_number, *phases[iter_number][]);
                    break;
                }
                iter_number++;
            }
        }
        //todo test
    };

    export
    struct AmplFreqOptExtractor {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        using InstFreqType = GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<InstFreqType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Central> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
            <double, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
                phase_computer;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
            <double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
                phase_computer_for_mode;

        INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
            decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                inst_freq_computer =
                    INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                        decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                            (integrator, derivator);

        INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq<double, decltype(integrator), decltype(derivator),
            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                inst_ampl_computer;
  
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, decltype(phase_computer)>
                inst_freq_computer_for_opt_filter = NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
                    NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, 
                        decltype(phase_computer)> (integrator, derivator, phase_computer);

        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, decltype(phase_computer_for_mode)>
                inst_freq_computer_for_mode_for_opt_filter = NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), 
                    decltype(derivator), NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, 
                        decltype(phase_computer_for_mode)> (integrator, derivator, phase_computer_for_mode);

        NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<double, decltype(integrator), decltype(derivator), decltype(inst_freq_computer_for_opt_filter),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                inst_ampl_computer_for_filtering = NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<double, decltype(integrator), 
                    decltype(derivator), decltype(inst_freq_computer_for_opt_filter),
                    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>(integrator, derivator, inst_freq_computer_for_opt_filter);
        
        NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter = NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);

        NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter2 = NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);
        
        NP_DSP::ONE_D::FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
            decltype(non_opt_filter2)> cascade_filter = NP_DSP::ONE_D::FILTERS::CascadeFilter
                <double, decltype(non_opt_filter), 
                    decltype(non_opt_filter2)>(non_opt_filter, non_opt_filter2);

        NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreq<double, decltype(integrator), decltype(derivator),
            decltype(cascade_filter), decltype(inst_ampl_computer_for_filtering)> req_filter = 
                NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreq<double, decltype(integrator), decltype(derivator),
                    decltype(cascade_filter), decltype(inst_ampl_computer_for_filtering)>(integrator, derivator, cascade_filter, inst_ampl_computer_for_filtering); 

        NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilterInstFreqDouble<double, decltype(req_filter), decltype(inst_freq_computer_for_opt_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)> filter2 =
                NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilter<double, decltype(req_filter), decltype(inst_freq_computer_for_opt_filter),
                    decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)>
                        (req_filter, inst_freq_computer_for_opt_filter, phase_computer,
                            inst_freq_computer_for_mode_for_opt_filter, phase_computer_for_mode);
        
        /*NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), decltype(inst_freq_computer_for_opt_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)> 
                filter2 = NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), 
                    decltype(inst_freq_computer_for_opt_filter), decltype(phase_computer), 
                        decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)> 
                            (cascade_filter, inst_freq_computer_for_opt_filter, phase_computer, 
                                inst_freq_computer_for_mode_for_opt_filter, phase_computer_for_mode);*/

        
        //FILTERS::CascadeFilter

        ~AmplFreqOptExtractor() {
            for (int i = 0; i < modes.size(); i++) {
                delete modes[i];
            }
            for (int i = 0; i < inst_freqs.size(); i++) {
                delete inst_freqs[i];
            }
            for (int i = 0; i < phases.size(); i++) {
                delete phases[i];
            }
            for (int i = 0; i < inst_ampls.size(); i++) {
                delete inst_ampls[i];
            }
        }

        template<Signal DataInT>
        void load(const DataInT& data_in) {
            data.base->vec->clear();
            for (auto i = 0; i < data_in.size(); i++) {
                data.base->vec->push_back(data_in[i]);
            }
            for (int i = 0; i < modes.size(); i++) {
                modes[i]->base->vec->clear();
            }
            for (int i = 0; i < inst_freqs.size(); i++) {
                inst_freqs[i]->base->vec->clear();
            }
            for (int i = 0; i < phases.size(); i++) {
                phases[i]->base->vec->clear();
            }
            for (int i = 0; i < inst_ampls.size(); i++) {
                inst_ampls[i]->base->vec->clear();
            }
        }

        bool stop(size_t iter_number) {
            IC(iter_number, (*phases[iter_number])[data.size()-1]);
            if ((*phases[iter_number])[data.size() - 1] < std::numbers::pi * 2.0) {
                return true;
            } else {
                return false;
            }
        }

        void compute() {
            inst_freq_computer.variability = 0.5;
            size_t iter_number = 0;
            for (;;) {
                //data.show(PlottingKind::Simple);
                if (iter_number + 1 > modes.size()) {
                    auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    modes.push_back(modes_new);
                    for (int i = 0; i < data.size(); i++) {
                        modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number + 1 > inst_freqs.size()) {
                    auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                    inst_freqs.push_back(inst_freqs_new);
                    for (int i = 0; i < data.size(); i++) {
                        inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
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
                        inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
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
                auto size_data = data.size();
                auto out_size = phases[iter_number]->size();
                phase_computer.compute(data, *phases[iter_number], nullptr);
                //
                phases[iter_number]->show(PlottingKind::Simple);

                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                filter2.compute(data, *modes[iter_number], inst_freqs[iter_number]);

                //
                modes[iter_number]->show(PlottingKind::Simple);
                for (int i = 0; i < data.size(); i++) {
                    (*modes[iter_number])[i] = data[i] - (*modes[iter_number])[i];
                }
                for (int i = 0; i < data.size(); i++) {
                    data[i] = data[i] - (*modes[iter_number])[i];
                }
                //
                modes[iter_number]->show(PlottingKind::Simple);

                DataType computer_buffer;
                for (int i = 0; i < data.size(); i++) {
                    computer_buffer.base->vec->push_back(0.0);
                }

                phase_computer_for_mode.compute(*modes[iter_number], *phases[iter_number], nullptr);
                //phases[iter_number]->show(PlottingKind::Simple);
                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                inst_ampl_computer.inst_freq_double = inst_freqs[iter_number];
                inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);

                //inst_ampls[iter_number]->show(PlottingKind::Simple);

                if (stop(iter_number)) {
                    iter_number++;

                    if (iter_number > modes.size() - 1) {
                        auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                        modes.push_back(modes_new);
                        for (int i = 0; i < data.size(); i++) {
                            modes[iter_number]->base->vec->push_back(0.0);
                        }
                    }
                    if (iter_number > inst_freqs.size() - 1) {
                        auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                        inst_freqs.push_back(inst_freqs_new);
                        for (int i = 0; i < data.size(); i++) {
                            inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                        }
                    }
                    if (iter_number > inst_ampls.size() - 1) {
                        auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                        inst_ampls.push_back(inst_ampls_new);
                        for (int i = 0; i < data.size(); i++) {
                            inst_ampls[iter_number]->base->vec->push_back(0.0);
                        }
                    }
                    if (iter_number > phases.size() - 1) {
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
                            inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
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

                    for (int i = 0; i < data.size(); i++) {
                        (*modes[iter_number])[i] = data[i];
                    }

                    phase_computer.compute(*modes[iter_number], *phases[iter_number], nullptr);
                    //phases[iter_number]->show(PlottingKind::Simple);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    //inst_freqs[iter_number]->show(PlottingKind::Simple);
                    inst_ampl_computer.inst_freq_double = inst_freqs[iter_number];
                    inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);
                    //inst_ampls[iter_number]->show(PlottingKind::Simple);
                    //IC(iter_number, *phases[iter_number][]);
                    break;
                }
                iter_number++;
            }
        }
        //todo test
    };

    export
    struct FreqAmplOptExtractor {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        using InstFreqType = GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<InstFreqType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Central> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
            <double, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
                phase_computer;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
            <double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
                phase_computer_for_mode;

        INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
            decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                inst_freq_computer =
                    INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                        decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                            (integrator, derivator);

        INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq<double, decltype(integrator), decltype(derivator),
            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                inst_ampl_computer;
  
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, decltype(phase_computer)>
                inst_freq_computer_for_opt_filter = NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
                    NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, 
                        decltype(phase_computer)> (integrator, derivator, phase_computer);

        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, decltype(phase_computer_for_mode)>
                inst_freq_computer_for_mode_for_opt_filter = NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), 
                    decltype(derivator), NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble, 
                        decltype(phase_computer_for_mode)> (integrator, derivator, phase_computer_for_mode);

        NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<double, decltype(integrator), decltype(derivator), decltype(inst_freq_computer_for_opt_filter),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                inst_ampl_computer_for_filtering = NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<double, decltype(integrator), 
                    decltype(derivator), decltype(inst_freq_computer_for_opt_filter),
                    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>(integrator, derivator, inst_freq_computer_for_opt_filter);
        
        NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter = NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);

        NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter2 = NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);
        
        NP_DSP::ONE_D::FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
            decltype(non_opt_filter2)> cascade_filter = NP_DSP::ONE_D::FILTERS::CascadeFilter
                <double, decltype(non_opt_filter), 
                    decltype(non_opt_filter2)>(non_opt_filter, non_opt_filter2);

        NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), decltype(inst_freq_computer_for_opt_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)> opt_filter =
                NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), decltype(inst_freq_computer_for_opt_filter),
                    decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)>
                        (cascade_filter, inst_freq_computer_for_opt_filter, phase_computer,
                            inst_freq_computer_for_mode_for_opt_filter, phase_computer_for_mode);

        NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreq<double, decltype(integrator), decltype(derivator),
            decltype(opt_filter), decltype(inst_ampl_computer_for_filtering)> filter2 = 
                NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreq<double, decltype(integrator), decltype(derivator),
                    decltype(opt_filter), decltype(inst_ampl_computer_for_filtering)>(integrator, derivator, opt_filter, inst_ampl_computer_for_filtering); 

        ~FreqAmplOptExtractor() {
            for (int i = 0; i < modes.size(); i++) {
                delete modes[i];
            }
            for (int i = 0; i < inst_freqs.size(); i++) {
                delete inst_freqs[i];
            }
            for (int i = 0; i < phases.size(); i++) {
                delete phases[i];
            }
            for (int i = 0; i < inst_ampls.size(); i++) {
                delete inst_ampls[i];
            }
        }

        template<Signal DataInT>
        void load(const DataInT& data_in) {
            data.base->vec->clear();
            for (auto i = 0; i < data_in.size(); i++) {
                data.base->vec->push_back(data_in[i]);
            }
            for (int i = 0; i < modes.size(); i++) {
                modes[i]->base->vec->clear();
            }
            for (int i = 0; i < inst_freqs.size(); i++) {
                inst_freqs[i]->base->vec->clear();
            }
            for (int i = 0; i < phases.size(); i++) {
                phases[i]->base->vec->clear();
            }
            for (int i = 0; i < inst_ampls.size(); i++) {
                inst_ampls[i]->base->vec->clear();
            }
        }

        bool stop(size_t iter_number) {
            IC(iter_number, (*phases[iter_number])[data.size()-1]);
            if ((*phases[iter_number])[data.size() - 1] < std::numbers::pi * 2.0) {
                return true;
            } else {
                return false;
            }
        }

        void compute() {
            inst_freq_computer.variability = 0.5;
            size_t iter_number = 0;
            for (;;) {
                //data.show(PlottingKind::Simple);
                if (iter_number + 1 > modes.size()) {
                    auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    modes.push_back(modes_new);
                    for (int i = 0; i < data.size(); i++) {
                        modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number + 1 > inst_freqs.size()) {
                    auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                    inst_freqs.push_back(inst_freqs_new);
                    for (int i = 0; i < data.size(); i++) {
                        inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
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
                        inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
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
                auto size_data = data.size();
                auto out_size = phases[iter_number]->size();
                phase_computer.compute(data, *phases[iter_number], nullptr);
                //
                phases[iter_number]->show(PlottingKind::Simple);

                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                filter2.compute(data, *modes[iter_number], inst_freqs[iter_number]);

                //
                modes[iter_number]->show(PlottingKind::Simple);
                for (int i = 0; i < data.size(); i++) {
                    (*modes[iter_number])[i] = data[i] - (*modes[iter_number])[i];
                }
                for (int i = 0; i < data.size(); i++) {
                    data[i] = data[i] - (*modes[iter_number])[i];
                }
                //
                modes[iter_number]->show(PlottingKind::Simple);

                DataType computer_buffer;
                for (int i = 0; i < data.size(); i++) {
                    computer_buffer.base->vec->push_back(0.0);
                }

                phase_computer_for_mode.compute(*modes[iter_number], *phases[iter_number], nullptr);
                //phases[iter_number]->show(PlottingKind::Simple);
                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                inst_ampl_computer.inst_freq_double = inst_freqs[iter_number];
                inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);

                //inst_ampls[iter_number]->show(PlottingKind::Simple);

                if (stop(iter_number)) {
                    iter_number++;

                    if (iter_number > modes.size() - 1) {
                        auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                        modes.push_back(modes_new);
                        for (int i = 0; i < data.size(); i++) {
                            modes[iter_number]->base->vec->push_back(0.0);
                        }
                    }
                    if (iter_number > inst_freqs.size() - 1) {
                        auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                        inst_freqs.push_back(inst_freqs_new);
                        for (int i = 0; i < data.size(); i++) {
                            inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                        }
                    }
                    if (iter_number > inst_ampls.size() - 1) {
                        auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                        inst_ampls.push_back(inst_ampls_new);
                        for (int i = 0; i < data.size(); i++) {
                            inst_ampls[iter_number]->base->vec->push_back(0.0);
                        }
                    }
                    if (iter_number > phases.size() - 1) {
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
                            inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
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

                    for (int i = 0; i < data.size(); i++) {
                        (*modes[iter_number])[i] = data[i];
                    }

                    phase_computer.compute(*modes[iter_number], *phases[iter_number], nullptr);
                    //phases[iter_number]->show(PlottingKind::Simple);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    //inst_freqs[iter_number]->show(PlottingKind::Simple);
                    inst_ampl_computer.inst_freq_double = inst_freqs[iter_number];
                    inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);
                    //inst_ampls[iter_number]->show(PlottingKind::Simple);
                    //IC(iter_number, *phases[iter_number][]);
                    break;
                }
                iter_number++;
            }
        }
        //todo test
    };
}
