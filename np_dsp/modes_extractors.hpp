#pragma once

#include <icecream.hpp>

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
    
    struct SimpleExtractorConfig {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<DataType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

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

        FILTERS::NonOptPeriodBasedFilter<double, FILTERS::FilteringType::AverageBased,
            decltype(integrator), FILTERS::InstFreqKind::Average>
                filter1;

        ~SimpleExtractorConfig() {
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
    };

    
    struct SimpleExtractorDoubleConfig {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        using InstFreqType = GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<InstFreqType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

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
            NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    filter1 = NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);

        

        ~SimpleExtractorDoubleConfig() {
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
    };

    
    struct AmplOptExtractorDoubleConfig {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        using InstFreqType = GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<InstFreqType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

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
            NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter = NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);

        NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreqDouble<double, decltype(integrator), decltype(derivator),
        decltype(non_opt_filter), decltype(inst_ampl_computer_for_filtering)> filter1 =
                NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreqDouble<double, decltype(integrator), decltype(derivator),
        decltype(non_opt_filter), decltype(inst_ampl_computer_for_filtering)> (integrator, derivator, non_opt_filter, inst_ampl_computer_for_filtering);

        //NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesDouble<double, decltype(integrator), decltype(derivator),
        //    decltype(non_opt_filter3), decltype(phase_computer), decltype(inst_ampl_computer_for_filtering), decltype(inst_freq_computer_for_mode_for_opt_filter)> 
        //        filter1 = 
        //        NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesDouble<double, decltype(integrator), decltype(derivator),
        //    decltype(non_opt_filter3), decltype(phase_computer), decltype(inst_ampl_computer_for_filtering), decltype(inst_freq_computer_for_mode_for_opt_filter)>
        //        (integrator, derivator, non_opt_filter3, phase_computer, inst_ampl_computer_for_filtering, inst_freq_computer_for_mode_for_opt_filter); 


        ~AmplOptExtractorDoubleConfig() {
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

        template<SignalBase DataInT>
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
    };

    
    struct FreqOptExtractorDoubleConfig {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        using InstFreqType = GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<InstFreqType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

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
            decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)> filter1 =
                NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), decltype(inst_freq_computer_for_opt_filter),
                    decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)>
                        (cascade_filter, inst_freq_computer_for_opt_filter, phase_computer,
                            inst_freq_computer_for_mode_for_opt_filter, phase_computer_for_mode);

        ~FreqOptExtractorDoubleConfig() {
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
        //todo test
    };
    
    
    struct FreqAmplOptExtractorDoubleConfig {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        using InstFreqType = GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<InstFreqType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

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

        NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesDouble<double, decltype(integrator), decltype(derivator),
            decltype(opt_filter), decltype(phase_computer), decltype(inst_ampl_computer_for_filtering), decltype(inst_freq_computer_for_mode_for_opt_filter)> filter1 = 
                NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesDouble<double, decltype(integrator), decltype(derivator),
            decltype(opt_filter), decltype(phase_computer), decltype(inst_ampl_computer_for_filtering), decltype(inst_freq_computer_for_mode_for_opt_filter)>
                (integrator, derivator, opt_filter, phase_computer, inst_ampl_computer_for_filtering, inst_freq_computer_for_mode_for_opt_filter);

        ~FreqAmplOptExtractorDoubleConfig() {
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
    };

    
    struct AmplOptExtractorConfig {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<DataType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

        PHASE_COMPUTERS::ExtremumsBasedNonOpt
            <double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
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

        INST_FREQ_COMPUTERS::PhaseBased
            <double, decltype(integrator), decltype(derivator), 
                INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                    decltype(phase_computer)> 
                        inst_freq_computer_for_filter =
        INST_FREQ_COMPUTERS::PhaseBased
            <double, decltype(integrator), decltype(derivator), 
                INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                    decltype(phase_computer)> 
                        (integrator, derivator, phase_computer);

        INST_FREQ_COMPUTERS::PhaseBased
            <double, decltype(integrator), decltype(derivator), 
                INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                    decltype(phase_computer_for_mode)> 
                        inst_freq_computer_for_filter_for_mode =
        INST_FREQ_COMPUTERS::PhaseBased
            <double, decltype(integrator), decltype(derivator), 
                INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                    decltype(phase_computer_for_mode)> 
                        (integrator, derivator, phase_computer_for_mode);

        INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<double, decltype(integrator), decltype(derivator), 
            decltype(inst_freq_computer_for_filter), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
            inst_ampl_computer_for_filter = 
        INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<double, decltype(integrator), decltype(derivator), 
            decltype(inst_freq_computer_for_filter), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                (integrator, derivator, inst_freq_computer_for_filter);

        FILTERS::NonOptPeriodBasedFilter<double, 
            FILTERS::FilteringType::DerivativeBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>
                    non_opt_filter = 
        FILTERS::NonOptPeriodBasedFilter<double, 
            FILTERS::FilteringType::DerivativeBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>(integrator);

        FILTERS::NonOptPeriodBasedFilter<double, 
            FILTERS::FilteringType::ValueBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>
                    non_opt_filter2 = 
        FILTERS::NonOptPeriodBasedFilter<double, 
            FILTERS::FilteringType::ValueBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>(integrator);

        FILTERS::NonOptPeriodBasedFilter<double, 
            FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>
                    non_opt_filter3 = 
        FILTERS::NonOptPeriodBasedFilter<double, 
            FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>(integrator);


        FILTERS::OptPeriodBasedFilter<double, decltype(non_opt_filter2), decltype(inst_freq_computer_for_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_filter_for_mode), decltype(phase_computer_for_mode)>
                opt_filter2 = 
        FILTERS::OptPeriodBasedFilter<double, decltype(non_opt_filter2), decltype(inst_freq_computer_for_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_filter_for_mode), decltype(phase_computer_for_mode)>
                (non_opt_filter2, inst_freq_computer_for_filter, phase_computer,
                    inst_freq_computer_for_filter_for_mode, phase_computer_for_mode);


        FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
            decltype(non_opt_filter2)> cascade_filter = 
        FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
            decltype(non_opt_filter2)>(non_opt_filter, non_opt_filter2);

        FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreq<double, decltype(integrator), decltype(derivator),
            decltype(cascade_filter), decltype(inst_ampl_computer_for_filter)> 
                filter = 
        FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreq<double, decltype(integrator), decltype(derivator),
            decltype(cascade_filter), decltype(inst_ampl_computer_for_filter)> 
                    (integrator, derivator, cascade_filter, inst_ampl_computer_for_filter);

                    FILTERS::NonOptPeriodBasedFilter<double, 
            FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>
                    filter1 = 
        FILTERS::NonOptPeriodBasedFilter<double, 
            FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>(integrator);
        
        /*FILTERS::OptPeriodBasedFilter<double, decltype(filter), decltype(inst_freq_computer_for_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_filter_for_mode), decltype(phase_computer_for_mode)>
                filter1 = 
        FILTERS::OptPeriodBasedFilter<double, decltype(filter), decltype(inst_freq_computer_for_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_filter_for_mode), decltype(phase_computer_for_mode)>
                (filter, inst_freq_computer_for_filter, phase_computer,
                    inst_freq_computer_for_filter_for_mode, phase_computer_for_mode);
        */

        //FILTERS::RecursiveFilterInstAmplChanges<double, decltype(integrator), decltype(derivator),
        //    decltype(non_opt_filter3), decltype(inst_freq_computer_for_filter), decltype(phase_computer), decltype(inst_ampl_computer_for_filter),
        //        decltype(inst_freq_computer_for_filter_for_mode), decltype(phase_computer_for_mode)> 
        //        filter1 = 
        //FILTERS::RecursiveFilterInstAmplChanges<double, decltype(integrator), decltype(derivator),
        //    decltype(non_opt_filter3), decltype(inst_freq_computer_for_filter), decltype(phase_computer), decltype(inst_ampl_computer_for_filter),
        //        decltype(inst_freq_computer_for_filter_for_mode), decltype(phase_computer_for_mode)>
        //            (integrator, derivator, non_opt_filter3, inst_freq_computer_for_filter, phase_computer, inst_ampl_computer_for_filter,
        //                inst_freq_computer_for_filter_for_mode, phase_computer_for_mode);
        
        
        
        
        ~AmplOptExtractorConfig() {
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
    };


    
    struct AmplOptExtractorDoubleAlternationFreqConfig {
        using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
        using InstFreqType = GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
        DataType data;
        std::vector<DataType *> modes;
        std::vector<InstFreqType *> inst_freqs;
        std::vector<DataType *> inst_ampls;
        std::vector<DataType *> phases;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

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

        NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreqDouble<double, decltype(integrator), decltype(derivator),
        decltype(cascade_filter), decltype(inst_ampl_computer_for_filtering)> filter1 =
                NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesFithConstInstFreqDouble<double, decltype(integrator), decltype(derivator),
        decltype(cascade_filter), decltype(inst_ampl_computer_for_filtering)> (integrator, derivator, cascade_filter, inst_ampl_computer_for_filtering);   
        
        /*NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesDouble<double, decltype(integrator), decltype(derivator),
            decltype(non_opt_filter2), decltype(phase_computer), decltype(inst_ampl_computer_for_filtering), decltype(inst_freq_computer_for_mode_for_opt_filter)> filter1 = 
                NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChangesDouble<double, decltype(integrator), decltype(derivator),
            decltype(non_opt_filter2), decltype(phase_computer), decltype(inst_ampl_computer_for_filtering), decltype(inst_freq_computer_for_mode_for_opt_filter)>
                (integrator, derivator, non_opt_filter2, phase_computer, inst_ampl_computer_for_filtering, inst_freq_computer_for_mode_for_opt_filter); */


        ~AmplOptExtractorDoubleAlternationFreqConfig() {
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
    };

     struct ByIterStopFunc{
        size_t max_iter_number = 50;
        template<typename Config>
        bool CheckStop(size_t iter_number, Config & config){
            if (iter_number > max_iter_number){
                return true;
            }
            return false;
        }
    };

     struct BySizeLogStopFunc{
        double log_muller = 1;

        template<typename Config>
        bool CheckStop(size_t iter_number, Config & config){
            if (std::log2(static_cast<double>(config.data.size())) * log_muller  < iter_number){
                return true;
            }
            return false;
        }
    };

     struct ByMaxPhaseStopFunc{
        double min_phase = 6.28;

        template<typename Config>
        bool CheckStop(size_t iter_number, Config & config){
            if ((*config.phases[config.phases.size() - 1])
                [(*config.phases[config.phases.size() - 1]).size() - 1] <= min_phase){
                return true;
            }
            return false;
        }
    };
    
    
    template<typename Config, typename Stop>
    void computeReq(Config & config, Stop & stop){
        size_t iter_number = 0;

        for (;;) {
            //data.show(PlottingKind::Simple);
            if (iter_number + 1 > config.modes.size()) {
                auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.modes.push_back(modes_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.modes[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.inst_freqs.size()) {
                auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.inst_freqs.push_back(inst_freqs_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_freqs[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.inst_ampls.size()) {
                auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.inst_ampls.push_back(inst_ampls_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.phases.size()) {
                auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.phases.push_back(phase_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.phases[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.modes[iter_number]->size() != config.data.size()) {
                config.modes[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.modes[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.inst_freqs[iter_number]->size() != config.data.size()) {
                config.inst_freqs[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_freqs[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.inst_ampls[iter_number]->size() != config.data.size()) {
                config.inst_ampls[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.phases[iter_number]->size() != config.data.size()) {
                config.phases[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.phases[iter_number]->base->vec->push_back(0.0);
                }
            }

            config.phase_computer.compute(config.data, *config.phases[iter_number], nullptr);
            //phases[iter_number]->show(PlottingKind::Simple);
            config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
            config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[iter_number]);
            //config.modes[iter_number]->show(PlottingKind::Simple);
            for (int i = 0; i < config.data.size(); i++) {
                (*config.modes[iter_number])[i] = config.data[i] - (*config.modes[iter_number])[i];
            }
            for (int i = 0; i < config.data.size(); i++) {
                config.data[i] = config.data[i] - (*config.modes[iter_number])[i];
            }

            if (stop.CheckStop(iter_number, config)) {
                typename Config::DataType computer_buffer;
                for (int i = 0; i < config.data.size(); i++) {
                    computer_buffer.base->vec->push_back(0.0);
                }
                iter_number++;
                if (iter_number > config.modes.size() - 1) {
                    auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.modes.push_back(modes_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.inst_freqs.size() - 1) {
                    auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.inst_freqs.push_back(inst_freqs_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_freqs[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.inst_ampls.size() - 1) {
                    auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.inst_ampls.push_back(inst_ampls_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.phases.size() - 1) {
                    auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.phases.push_back(phase_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.modes[iter_number]->size() != config.data.size()) {
                    config.modes[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.inst_freqs[iter_number]->size() != config.data.size()) {
                    config.inst_freqs[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_freqs[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.inst_ampls[iter_number]->size() != config.data.size()) {
                    config.inst_ampls[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.phases[iter_number]->size() != config.data.size()) {
                    config.phases[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                for (int i = 0; i < config.data.size(); i++) {
                    (*config.modes[iter_number])[i] = config.data[i];
                }
                for(int i = 0; i < iter_number; i++){
                    config.phase_computer_for_mode.compute(*config.modes[i], *config.phases[i], nullptr);
                    config.inst_freq_computer.compute(*config.phases[i], *config.inst_freqs[i], nullptr);
                    config.inst_ampl_computer.inst_freq = config.inst_freqs[i];
                    config.inst_ampl_computer.compute(*config.modes[i], *config.inst_ampls[i], &computer_buffer);
                }
                break;
            }
            iter_number++;
        }
    }

    
    template<typename Config, typename Stop>
    void computeReqDouble(Config & config, Stop & stop){
        size_t iter_number = 0;
        for (;;) {
            //data.show(PlottingKind::Simple);
            if (iter_number + 1 > config.modes.size()) {
                auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.modes.push_back(modes_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.modes[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.inst_freqs.size()) {
                auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                config.inst_freqs.push_back(inst_freqs_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                }
            }
            if (iter_number + 1 > config.inst_ampls.size()) {
                auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.inst_ampls.push_back(inst_ampls_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.phases.size()) {
                auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.phases.push_back(phase_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.phases[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.modes[iter_number]->size() != config.data.size()) {
                config.modes[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.modes[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.inst_freqs[iter_number]->size() != config.data.size()) {
                config.inst_freqs[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                }
            }
            if (config.inst_ampls[iter_number]->size() != config.data.size()) {
                config.inst_ampls[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.phases[iter_number]->size() != config.data.size()) {
                config.phases[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.phases[iter_number]->base->vec->push_back(0.0);
                }
            }
            /*for (int i = 0; i + 1 < iter_number; i++){
                IC(iter_number, i);
                config.inst_freq_computer.compute(*config.phases[i], *config.inst_freqs[i], nullptr);
                config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[i]);
                for(int j = 0; j < config.data.size(); j++){
                    (*config.modes[i])[j] = (*config.modes[i])[j] + (config.data[j] - (*config.modes[iter_number])[j]);
                    config.data[j] = (*config.modes[iter_number])[j];
                }
                
                config.inst_freq_computer.compute(*config.phases[i], *config.inst_freqs[i], nullptr);
                config.filter1.compute(*config.modes[i], *config.modes[iter_number], config.inst_freqs[i]);
                for(int j = 0; j < config.data.size(); j++){
                    (*config.modes[i])[j] = (*config.modes[i])[j] - (*config.modes[iter_number])[j];
                    (*config.modes[i+1])[j] = (*config.modes[i+1])[j] + (*config.modes[iter_number])[j];
                }
                
                if(i!=0){
                    config.inst_freq_computer.compute(*config.phases[i - 1], *config.inst_freqs[i - 1], nullptr);
                    config.filter1.compute(*config.modes[i], *config.modes[iter_number], config.inst_freqs[i -1]);
                    for(int j = 0; j < config.data.size(); j++){
                        (*config.modes[i-1])[j] = (*config.modes[i-1])[j] + ((*config.modes[i])[j] - (*config.modes[iter_number])[j]);
                        (*config.modes[i])[j] = (*config.modes[iter_number])[j];
                    }
                }
            }*/
            
            config.phase_computer.compute(config.data, *config.phases[iter_number], nullptr);
            //phases[iter_number]->show(PlottingKind::Simple);
            config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
            config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[iter_number]);
            //config.modes[iter_number]->show(PlottingKind::Simple);
            for (int i = 0; i < config.data.size(); i++) {
                (*config.modes[iter_number])[i] = config.data[i] - (*config.modes[iter_number])[i];
            }
            for (int i = 0; i < config.data.size(); i++) {
                config.data[i] = config.data[i] - (*config.modes[iter_number])[i];
            }
            //config.modes[iter_number]->show(PlottingKind::Simple);
            
            //phase_computer_for_mode.compute(*modes[iter_number], *phases[iter_number], nullptr);
            //phases[iter_number]->show(PlottingKind::Simple);
            //inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
            //inst_ampl_computer.inst_freq = inst_freqs[iter_number];
            //inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);
            //inst_ampls[iter_number]->show(PlottingKind::Simple);
            if (stop.CheckStop(iter_number, config)) {
                typename Config::DataType computer_buffer;
                for (int i = 0; i < config.data.size(); i++) {
                    computer_buffer.base->vec->push_back(0.0);
                }
                iter_number++;
                if (iter_number > config.modes.size() - 1) {
                    auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.modes.push_back(modes_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.inst_freqs.size() - 1) {
                    auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                    config.inst_freqs.push_back(inst_freqs_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                    }
                }
                if (iter_number > config.inst_ampls.size() - 1) {
                    auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.inst_ampls.push_back(inst_ampls_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.phases.size() - 1) {
                    auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.phases.push_back(phase_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.modes[iter_number]->size() != config.data.size()) {
                    config.modes[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.inst_freqs[iter_number]->size() != config.data.size()) {
                    config.inst_freqs[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                    }
                }
                if (config.inst_ampls[iter_number]->size() != config.data.size()) {
                    config.inst_ampls[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.phases[iter_number]->size() != config.data.size()) {
                    config.phases[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                for (int i = 0; i < config.data.size(); i++) {
                    (*config.modes[iter_number])[i] = config.data[i];
                }
                for(int i = 0; i < iter_number; i++){
                    config.phase_computer_for_mode.compute(*config.modes[i], *config.phases[i], nullptr);
                    config.inst_freq_computer.compute(*config.phases[i], *config.inst_freqs[i], nullptr);
                    config.inst_ampl_computer.inst_freq_double = config.inst_freqs[i];
                    config.inst_ampl_computer.compute(*config.modes[i], *config.inst_ampls[i], &computer_buffer);
                }
                break;
            }
            iter_number++;
        }
    }

    
    template<typename Config, typename Stop>
    void computeReqCascade(Config & config, Stop & stop){
        size_t iter_number = 0;

        for (;;) {
            //data.show(PlottingKind::Simple);
            if (iter_number + 1 > config.modes.size()) {
                auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.modes.push_back(modes_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.modes[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.inst_freqs.size()) {
                auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.inst_freqs.push_back(inst_freqs_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_freqs[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.inst_ampls.size()) {
                auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.inst_ampls.push_back(inst_ampls_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.phases.size()) {
                auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.phases.push_back(phase_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.phases[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.modes[iter_number]->size() != config.data.size()) {
                config.modes[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.modes[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.inst_freqs[iter_number]->size() != config.data.size()) {
                config.inst_freqs[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_freqs[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.inst_ampls[iter_number]->size() != config.data.size()) {
                config.inst_ampls[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.phases[iter_number]->size() != config.data.size()) {
                config.phases[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.phases[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number % 2 == 0){
                config.phase_computer.compute(config.data, *config.phases[iter_number], nullptr);
                //phases[iter_number]->show(PlottingKind::Simple);
                config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
                config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[iter_number]);
                //config.modes[iter_number]->show(PlottingKind::Simple);
                for (int i = 0; i < config.data.size(); i++) {
                    (*config.modes[iter_number])[i] = config.data[i] - (*config.modes[iter_number])[i];
                }
                for (int i = 0; i < config.data.size(); i++) {
                    config.data[i] = config.data[i] - (*config.modes[iter_number])[i];
                }
            }
            else{
                config.phase_computer.compute(*config.modes[iter_number-1], *config.phases[iter_number], nullptr);
                //phases[iter_number]->show(PlottingKind::Simple);
                config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
                config.filter1.compute(*config.modes[iter_number-1], *config.modes[iter_number], config.inst_freqs[iter_number]);
                //config.modes[iter_number]->show(PlottingKind::Simple);
                for (int i = 0; i < config.data.size(); i++) {
                    auto bf = (*config.modes[iter_number - 1])[i];
                    (*config.modes[iter_number - 1])[i] = (*config.modes[iter_number])[i];
                    (*config.modes[iter_number])[i] = bf - (*config.modes[iter_number - 1])[i];
                }
            }

            

            if (stop.CheckStop(iter_number, config)) {
                typename Config::DataType computer_buffer;
                for (int i = 0; i < config.data.size(); i++) {
                    computer_buffer.base->vec->push_back(0.0);
                }
                iter_number++;
                if (iter_number > config.modes.size() - 1) {
                    auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.modes.push_back(modes_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.inst_freqs.size() - 1) {
                    auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.inst_freqs.push_back(inst_freqs_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_freqs[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.inst_ampls.size() - 1) {
                    auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.inst_ampls.push_back(inst_ampls_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.phases.size() - 1) {
                    auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.phases.push_back(phase_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.modes[iter_number]->size() != config.data.size()) {
                    config.modes[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.inst_freqs[iter_number]->size() != config.data.size()) {
                    config.inst_freqs[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_freqs[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.inst_ampls[iter_number]->size() != config.data.size()) {
                    config.inst_ampls[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.phases[iter_number]->size() != config.data.size()) {
                    config.phases[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                for (int i = 0; i < config.data.size(); i++) {
                    (*config.modes[iter_number])[i] = config.data[i];
                }
                for(int i = 0; i < iter_number; i++){
                    config.phase_computer_for_mode.compute(*config.modes[i], *config.phases[i], nullptr);
                    config.inst_freq_computer.compute(*config.phases[i], *config.inst_freqs[i], nullptr);
                    config.inst_ampl_computer.inst_freq = config.inst_freqs[i];
                    config.inst_ampl_computer.compute(*config.modes[i], *config.inst_ampls[i], &computer_buffer);
                }
                break;
            }
            iter_number++;
        }
    }

    
    template<typename Config, typename Stop>
    void computeReqDoubleCascade(Config & config, Stop & stop){
        size_t iter_number = 0;
        for (;;) {
            //data.show(PlottingKind::Simple);
            if (iter_number + 1 > config.modes.size()) {
                auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.modes.push_back(modes_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.modes[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.inst_freqs.size()) {
                auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                config.inst_freqs.push_back(inst_freqs_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                }
            }
            if (iter_number + 1 > config.inst_ampls.size()) {
                auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.inst_ampls.push_back(inst_ampls_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (iter_number + 1 > config.phases.size()) {
                auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                config.phases.push_back(phase_new);
                for (int i = 0; i < config.data.size(); i++) {
                    config.phases[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.modes[iter_number]->size() != config.data.size()) {
                config.modes[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.modes[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.inst_freqs[iter_number]->size() != config.data.size()) {
                config.inst_freqs[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                }
            }
            if (config.inst_ampls[iter_number]->size() != config.data.size()) {
                config.inst_ampls[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                }
            }
            if (config.phases[iter_number]->size() != config.data.size()) {
                config.phases[iter_number]->base->vec->clear();
                for (int i = 0; i < config.data.size(); i++) {
                    config.phases[iter_number]->base->vec->push_back(0.0);
                }
            }
            /*for (int i = 0; i + 1 < iter_number; i++){
                IC(iter_number, i);
                config.inst_freq_computer.compute(*config.phases[i], *config.inst_freqs[i], nullptr);
                config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[i]);
                for(int j = 0; j < config.data.size(); j++){
                    (*config.modes[i])[j] = (*config.modes[i])[j] + (config.data[j] - (*config.modes[iter_number])[j]);
                    config.data[j] = (*config.modes[iter_number])[j];
                }
                
                config.inst_freq_computer.compute(*config.phases[i], *config.inst_freqs[i], nullptr);
                config.filter1.compute(*config.modes[i], *config.modes[iter_number], config.inst_freqs[i]);
                for(int j = 0; j < config.data.size(); j++){
                    (*config.modes[i])[j] = (*config.modes[i])[j] - (*config.modes[iter_number])[j];
                    (*config.modes[i+1])[j] = (*config.modes[i+1])[j] + (*config.modes[iter_number])[j];
                }
                
                if(i!=0){
                    config.inst_freq_computer.compute(*config.phases[i - 1], *config.inst_freqs[i - 1], nullptr);
                    config.filter1.compute(*config.modes[i], *config.modes[iter_number], config.inst_freqs[i -1]);
                    for(int j = 0; j < config.data.size(); j++){
                        (*config.modes[i-1])[j] = (*config.modes[i-1])[j] + ((*config.modes[i])[j] - (*config.modes[iter_number])[j]);
                        (*config.modes[i])[j] = (*config.modes[iter_number])[j];
                    }
                }
            }*/

            if (iter_number % 2 == 0){
                config.phase_computer.compute(config.data, *config.phases[iter_number], nullptr);
                //phases[iter_number]->show(PlottingKind::Simple);
                config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
                config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[iter_number]);
                //config.modes[iter_number]->show(PlottingKind::Simple);
                for (int i = 0; i < config.data.size(); i++) {
                    (*config.modes[iter_number])[i] = config.data[i] - (*config.modes[iter_number])[i];
                }
                for (int i = 0; i < config.data.size(); i++) {
                    config.data[i] = config.data[i] - (*config.modes[iter_number])[i];
                }
            }
            else{
                config.phase_computer.compute(*config.modes[iter_number-1], *config.phases[iter_number], nullptr);
                //phases[iter_number]->show(PlottingKind::Simple);
                config.inst_freq_computer.compute(*config.phases[iter_number-1], *config.inst_freqs[iter_number], nullptr);
                config.filter1.compute(*config.modes[iter_number-1], *config.modes[iter_number], config.inst_freqs[iter_number]);
                //config.modes[iter_number]->show(PlottingKind::Simple);
                for (int i = 0; i < config.data.size(); i++) {
                    auto bf = (*config.modes[iter_number - 1])[i];
                    (*config.modes[iter_number - 1])[i] = (*config.modes[iter_number])[i];
                    (*config.modes[iter_number])[i] = bf - (*config.modes[iter_number - 1])[i];
                }
            }
            //config.modes[iter_number]->show(PlottingKind::Simple);
            
            //phase_computer_for_mode.compute(*modes[iter_number], *phases[iter_number], nullptr);
            //phases[iter_number]->show(PlottingKind::Simple);
            //inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
            //inst_ampl_computer.inst_freq = inst_freqs[iter_number];
            //inst_ampl_computer.compute(*modes[iter_number], *inst_ampls[iter_number], &computer_buffer);
            //inst_ampls[iter_number]->show(PlottingKind::Simple);
            if (stop.CheckStop(iter_number, config)) {
                typename Config::DataType computer_buffer;
                for (int i = 0; i < config.data.size(); i++) {
                    computer_buffer.base->vec->push_back(0.0);
                }
                iter_number++;
                if (iter_number > config.modes.size() - 1) {
                    auto* modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.modes.push_back(modes_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.inst_freqs.size() - 1) {
                    auto* inst_freqs_new = new GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true>;
                    config.inst_freqs.push_back(inst_freqs_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                    }
                }
                if (iter_number > config.inst_ampls.size() - 1) {
                    auto* inst_ampls_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.inst_ampls.push_back(inst_ampls_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (iter_number > config.phases.size() - 1) {
                    auto* phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
                    config.phases.push_back(phase_new);
                    for (int i = 0; i < config.data.size(); i++) {
                        config.phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.modes[iter_number]->size() != config.data.size()) {
                    config.modes[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.modes[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.inst_freqs[iter_number]->size() != config.data.size()) {
                    config.inst_freqs[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_freqs[iter_number]->base->vec->push_back({0.0, 0.0});
                    }
                }
                if (config.inst_ampls[iter_number]->size() != config.data.size()) {
                    config.inst_ampls[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.inst_ampls[iter_number]->base->vec->push_back(0.0);
                    }
                }
                if (config.phases[iter_number]->size() != config.data.size()) {
                    config.phases[iter_number]->base->vec->clear();
                    for (int i = 0; i < config.data.size(); i++) {
                        config.phases[iter_number]->base->vec->push_back(0.0);
                    }
                }
                for (int i = 0; i < config.data.size(); i++) {
                    (*config.modes[iter_number])[i] = config.data[i];
                }
                for(int i = 0; i < iter_number; i++){
                    config.phase_computer_for_mode.compute(*config.modes[i], *config.phases[i], nullptr);
                    config.inst_freq_computer.compute(*config.phases[i], *config.inst_freqs[i], nullptr);
                    config.inst_ampl_computer.inst_freq_double = config.inst_freqs[i];
                    config.inst_ampl_computer.compute(*config.modes[i], *config.inst_ampls[i], &computer_buffer);
                }
                break;
            }
            iter_number++;
        }
    }

}
