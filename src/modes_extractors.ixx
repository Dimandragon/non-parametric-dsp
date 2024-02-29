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
    template<InstFreqComputer InstFreqComputerT>
    struct LinearExtractor {
        //using InstFreqComputerT
    };

    export
    struct MainExtractor {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<DataType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<DataType, DataType, INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DataType, DataType, DERIVATORS::FinniteDifferenceType::Central> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
        phase_computer;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
        phase_computer_for_mode;

        INST_FREQ_COMPUTERS::ComputedOnPhase<DataType, DataType, decltype(integrator),
            decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
        inst_freq_computer =
                INST_FREQ_COMPUTERS::ComputedOnPhase<DataType, DataType, decltype(integrator),
                    decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                (integrator, derivator);

        INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq<DataType, DataType, DataType,
            DataType, decltype(integrator), decltype(derivator),
            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
        inst_ampl_computer;


        /*FILTERS::OptPeriodBasedFilter<DataType, DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::DerArctg,
            FILTERS::FilteringType::DerivativeBased, decltype(integrator), decltype(derivator), FILTERS::InstFreqComputerKind::phase_based_time_average,
                FILTERS::PhaseComputingKind::extremums_based_non_opt>
                    filter = FILTERS::OptPeriodBasedFilter<DataType, DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::DerArctg, FILTERS::FilteringType::DerivativeBased,
                        decltype(integrator), decltype(derivator), FILTERS::InstFreqComputerKind::phase_based_time_average,
                            FILTERS::PhaseComputingKind::extremums_based_non_opt>(integrator, derivator);*/

        FILTERS::NonOptPeriodBasedFilter<DataType, DataType,
            DataType, FILTERS::FilteringType::Median,
            decltype(integrator), FILTERS::InstFreqKind::Average>
        filter1;

        FILTERS::NonOptPeriodBasedFilter<DataType, DataType,
            DataType, FILTERS::FilteringType::DerivativeBased,
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
                inst_ampls[iter_number]->kind = SignalKind::Stohastic;
                modes[iter_number]->kind = SignalKind::Stohastic;
                inst_freqs[iter_number]->kind = SignalKind::Stohastic;
                data.kind = SignalKind::Stohastic;
                phases[iter_number]->kind = SignalKind::Monotone;
                GENERAL::Nil nil;

                phase_computer.compute(data, *phases[iter_number], nil);
                //phases[iter_number]->show(PlottingKind::Simple);

                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nil);
                filter2.compute(data, *modes[iter_number], *inst_freqs[iter_number]);


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
                phase_computer_for_mode.compute(*modes[iter_number], *phases[iter_number], nil);
                //phases[iter_number]->show(PlottingKind::Simple);
                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nil);
                inst_ampl_computer.inst_freq = inst_freqs[iter_number];
                inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], computer_buffer);

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
                    inst_ampls[iter_number]->kind = SignalKind::Stohastic;
                    modes[iter_number]->kind = SignalKind::Stohastic;
                    inst_freqs[iter_number]->kind = SignalKind::Stohastic;
                    data.kind = SignalKind::Stohastic;
                    phases[iter_number]->kind = SignalKind::Monotone;


                    phase_computer.compute(*modes[iter_number], *phases[iter_number], nil);
                    //phases[iter_number]->show(PlottingKind::Simple);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nil);
                    //inst_freqs[iter_number]->show(PlottingKind::Simple);
                    inst_ampl_computer.inst_freq = inst_freqs[iter_number];
                    inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], computer_buffer);
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

        INTEGRATORS::Riman<DataType, DataType, INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DataType, DataType, DERIVATORS::FinniteDifferenceType::Central> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
        phase_computer;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
        phase_computer_for_mode;

        INST_FREQ_COMPUTERS::ComputedOnPhase<DataType, InstFreqType, decltype(integrator),
            decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
        inst_freq_computer =
                INST_FREQ_COMPUTERS::ComputedOnPhase<DataType, InstFreqType, decltype(integrator),
                    decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
                (integrator, derivator);

        INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq<DataType, DataType, DataType,
            InstFreqType, decltype(integrator), decltype(derivator),
            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveDouble>
        inst_ampl_computer;


        /*FILTERS::OptPeriodBasedFilter<DataType, DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::DerArctg,
            FILTERS::FilteringType::DerivativeBased, decltype(integrator), decltype(derivator), FILTERS::InstFreqComputerKind::phase_based_time_average,
                FILTERS::PhaseComputingKind::extremums_based_non_opt>
                    filter = FILTERS::OptPeriodBasedFilter<DataType, DataType, DataType, PHASE_COMPUTERS::ExtremumsKind::DerArctg, FILTERS::FilteringType::DerivativeBased,
                        decltype(integrator), decltype(derivator), FILTERS::InstFreqComputerKind::phase_based_time_average,
                            FILTERS::PhaseComputingKind::extremums_based_non_opt>(integrator, derivator);*/

        FILTERS::NonOptPeriodBasedFilter<DataType, DataType,
            InstFreqType, FILTERS::FilteringType::Median,
            decltype(integrator), FILTERS::InstFreqKind::Double>
        filter1;

        FILTERS::NonOptPeriodBasedFilter<DataType, DataType,
            InstFreqType, FILTERS::FilteringType::DerivativeBased,
            decltype(integrator), FILTERS::InstFreqKind::Double>
        filter2;

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
                GENERAL::Nil nil;
                data.kind = SignalKind::Stohastic;
                phases[iter_number]->kind = SignalKind::Monotone;
                phase_computer.compute(data, *phases[iter_number], nil);
                //phases[iter_number]->show(PlottingKind::Simple);

                inst_freqs[iter_number]->kind = SignalKind::Stohastic;
                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nil);
                //inst_freqs[iter_number]->show(PlottingKind::Simple);

                //if (iter_number%2 == 0) {
                //filter1.compute(data, *modes[iter_number], *inst_freqs[iter_number]);
                //}
                //else {
                modes[iter_number]->kind = SignalKind::Stohastic;
                filter2.compute(data, *modes[iter_number], *inst_freqs[iter_number]);
                //}


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

                phase_computer_for_mode.compute(*modes[iter_number], *phases[iter_number], nil);
                //phases[iter_number]->show(PlottingKind::Simple);
                inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nil);
                inst_ampl_computer.inst_freq = inst_freqs[iter_number];
                inst_ampls[iter_number]->kind = SignalKind::Stohastic;
                inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], computer_buffer);

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
                    inst_ampls[iter_number]->kind = SignalKind::Stohastic;
                    modes[iter_number]->kind = SignalKind::Stohastic;
                    inst_freqs[iter_number]->kind = SignalKind::Stohastic;
                    data.kind = SignalKind::Stohastic;
                    phases[iter_number]->kind = SignalKind::Monotone;

                    phase_computer.compute(*modes[iter_number], *phases[iter_number], nil);
                    //phases[iter_number]->show(PlottingKind::Simple);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nil);
                    //inst_freqs[iter_number]->show(PlottingKind::Simple);
                    inst_ampl_computer.inst_freq = inst_freqs[iter_number];
                    inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], computer_buffer);
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
