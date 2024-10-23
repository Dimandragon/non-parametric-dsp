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

        FILTERS::InstFreqBased<double, FILTERS::FilteringType::AverageBased,
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
        
        NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    filter1 = NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
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

        NP_DSP::ONE_D::INST_AMPL_COMPUTERS::HilbertTransformBased
            <NP_DSP::ONE_D::UTILITY_MATH::HTKind::Mull> inst_ampl_computer;

        NP_DSP::ONE_D::INST_AMPL_COMPUTERS::HilbertTransformBased
            <NP_DSP::ONE_D::UTILITY_MATH::HTKind::Mull> inst_ampl_computer_for_filtering;
  
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

        NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter = NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);

        NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::RecursiveFilterInstAmplChangesWithConstInstFreqDouble<double, decltype(integrator), decltype(derivator),
        decltype(non_opt_filter), decltype(inst_ampl_computer_for_filtering)> filter1 =
                NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::RecursiveFilterInstAmplChangesWithConstInstFreqDouble<double, decltype(integrator), decltype(derivator),
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
        
        NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter = NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);

        NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter2 = NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);
        
        NP_DSP::ONE_D::FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
            decltype(non_opt_filter2)> cascade_filter = NP_DSP::ONE_D::FILTERS::CascadeFilter
                <double, decltype(non_opt_filter), 
                    decltype(non_opt_filter2)>(non_opt_filter, non_opt_filter2);

        NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), decltype(inst_freq_computer_for_opt_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)> filter1 =
                NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), decltype(inst_freq_computer_for_opt_filter),
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
        
        NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter = NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);

        NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter2 = NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);
        
        NP_DSP::ONE_D::FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
            decltype(non_opt_filter2)> cascade_filter = NP_DSP::ONE_D::FILTERS::CascadeFilter
                <double, decltype(non_opt_filter), 
                    decltype(non_opt_filter2)>(non_opt_filter, non_opt_filter2);

        NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), decltype(inst_freq_computer_for_opt_filter),
            decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)> opt_filter =
                NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::OptPeriodBasedFilterInstFreqDouble<double, decltype(cascade_filter), decltype(inst_freq_computer_for_opt_filter),
                    decltype(phase_computer), decltype(inst_freq_computer_for_mode_for_opt_filter), decltype(phase_computer_for_mode)>
                        (cascade_filter, inst_freq_computer_for_opt_filter, phase_computer,
                            inst_freq_computer_for_mode_for_opt_filter, phase_computer_for_mode);

        NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::RecursiveFilterInstAmplChangesDouble<double, decltype(integrator), decltype(derivator),
            decltype(opt_filter), decltype(phase_computer), decltype(inst_ampl_computer_for_filtering), decltype(inst_freq_computer_for_mode_for_opt_filter)> filter1 = 
                NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::RecursiveFilterInstAmplChangesDouble<double, decltype(integrator), decltype(derivator),
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

        /*INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq<double, decltype(integrator), decltype(derivator),
            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_ampl_computer;*/
        NP_DSP::ONE_D::INST_AMPL_COMPUTERS::HilbertTransformBased
            <NP_DSP::ONE_D::UTILITY_MATH::HTKind::Mull> inst_ampl_computer;

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

        NP_DSP::ONE_D::INST_AMPL_COMPUTERS::HilbertTransformBased
            <NP_DSP::ONE_D::UTILITY_MATH::HTKind::Mull>
                inst_ampl_computer_for_filter;

        FILTERS::InstFreqBased<double, 
            FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>
                    non_opt_filter = 
        FILTERS::InstFreqBased<double, 
            FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>(integrator);

        FILTERS::EXPERIMENTAL::RecursiveFilterInstAmplChangesWithConstInstFreq<double, decltype(integrator), decltype(derivator),
            decltype(non_opt_filter), decltype(inst_ampl_computer_for_filter)> 
                filter2 = 
        FILTERS::EXPERIMENTAL::RecursiveFilterInstAmplChangesWithConstInstFreq<double, decltype(integrator), decltype(derivator),
            decltype(non_opt_filter), decltype(inst_ampl_computer_for_filter)> 
                (integrator, derivator, non_opt_filter, inst_ampl_computer_for_filter);      
        
        
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

        NP_DSP::ONE_D::INST_AMPL_COMPUTERS::HilbertTransformBased
            <NP_DSP::ONE_D::UTILITY_MATH::HTKind::Mull> inst_ampl_computer;

        NP_DSP::ONE_D::INST_AMPL_COMPUTERS::HilbertTransformBased
            <NP_DSP::ONE_D::UTILITY_MATH::HTKind::Mull> inst_ampl_computer_for_filtering;

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

        NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
            NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
                decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                    non_opt_filter = NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
                        NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
                            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Double>
                                (integrator);

        NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::RecursiveFilterInstAmplChangesWithConstInstFreqDouble<double, decltype(integrator), decltype(derivator),
        decltype(non_opt_filter), decltype(inst_ampl_computer_for_filtering)> filter1 =
                NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::RecursiveFilterInstAmplChangesWithConstInstFreqDouble<double, decltype(integrator), decltype(derivator),
        decltype(non_opt_filter), decltype(inst_ampl_computer_for_filtering)> (integrator, derivator, non_opt_filter, inst_ampl_computer_for_filtering);   
        
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
            config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
            config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[iter_number]);
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
            
            
            config.phase_computer.compute(config.data, *config.phases[iter_number], nullptr);
            config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
            config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[iter_number]);
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
                    if constexpr(decltype(config.inst_ampl_computer)::is_used_external_inst_freq){
                        config.inst_ampl_computer.inst_freq_double = config.inst_freqs[i];
                    }
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
                config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
                config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[iter_number]);
                for (int i = 0; i < config.data.size(); i++) {
                    (*config.modes[iter_number])[i] = config.data[i] - (*config.modes[iter_number])[i];
                }
                for (int i = 0; i < config.data.size(); i++) {
                    config.data[i] = config.data[i] - (*config.modes[iter_number])[i];
                }
            }
            else{
                config.phase_computer.compute(*config.modes[iter_number-1], *config.phases[iter_number], nullptr);
                config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
                config.filter1.compute(*config.modes[iter_number-1], *config.modes[iter_number], config.inst_freqs[iter_number]);
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
            if (iter_number % 2 == 0){
                config.phase_computer.compute(config.data, *config.phases[iter_number], nullptr);
                config.inst_freq_computer.compute(*config.phases[iter_number], *config.inst_freqs[iter_number], nullptr);
                config.filter1.compute(config.data, *config.modes[iter_number], config.inst_freqs[iter_number]);
                for (int i = 0; i < config.data.size(); i++) {
                    (*config.modes[iter_number])[i] = config.data[i] - (*config.modes[iter_number])[i];
                }
                for (int i = 0; i < config.data.size(); i++) {
                    config.data[i] = config.data[i] - (*config.modes[iter_number])[i];
                }
            }
            else{
                config.phase_computer.compute(*config.modes[iter_number-1], *config.phases[iter_number], nullptr);
                config.inst_freq_computer.compute(*config.phases[iter_number-1], *config.inst_freqs[iter_number], nullptr);
                config.filter1.compute(*config.modes[iter_number-1], *config.modes[iter_number], config.inst_freqs[iter_number]);
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
                    if constexpr(decltype(config.inst_ampl_computer)::is_used_external_inst_freq){
                        config.inst_ampl_computer.inst_freq_double = config.inst_freqs[i];
                    }
                    config.inst_ampl_computer.compute(*config.modes[i], *config.inst_ampls[i], &computer_buffer);
                }
                break;
            }
            iter_number++;
        }
    }

    struct InstFreqAndAmplNormSincExtractor
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

        FILTERS::InstFreqBased<double,
                FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>
        non_opt_filter = FILTERS::InstFreqBased<double,
                FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>(integrator);

        FILTERS::MonoFreqFilters<double, FILTERS::MonoInstFreqFilteringType::SincPaddedFIR> filter;

        template<Signal DataT>
        void compute(const DataT & data_in){
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

            while(true){
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

                phase_computer_simple.compute(data, *phases[iter_number], nullptr);

                if((*phases[iter_number])[data.size() - 1] > 6.28){
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    double base_inst_freq = INST_FREQ_COMPUTERS::instFreqNorm(data, data_buffer, *inst_freqs[iter_number], freq_conv, freq_conv_image);

                    phase_computer_simple.compute(data_buffer, *phases[iter_number], nullptr);

                    base_inst_freq = 1.0 /
                            (static_cast<double>(data.size()) /
                                ((*phases[iter_number])[data.size() - 1] / 2.0 / std::numbers::pi));

                    for(int i = 0; i < data.size(); i++){
                        data[i] = data_buffer[i];
                    }
                    filter.is_low_pass = false;
                    filter.freq = base_inst_freq;
                    filter.compute(data, *modes[iter_number], *phases[iter_number]); // phases is unused


                    filter.is_low_pass = true;
                    filter.compute(data, data_buffer, *phases[iter_number]);

                    auto size_expr = [&](){
                        return data.size();
                    };

                    auto val_expr = [&](size_t idx){
                        return base_inst_freq;
                    };

                    ExpressionWrapper<double, size_t, decltype(val_expr),
                            GENERAL::Nil, decltype(size_expr), false> inst_freq_buffer(val_expr, size_expr);

                    FILTERS::EXPERIMENTAL::InstAmplNormalizatorNaive<double, decltype(inst_ampl_computer),
                            decltype(inst_freq_buffer), decltype(non_opt_filter)>
                            naive_normalizer = FILTERS::EXPERIMENTAL::InstAmplNormalizatorNaive<double, decltype(inst_ampl_computer),
                            decltype(inst_freq_buffer), decltype(non_opt_filter)>
                            (inst_ampl_computer, non_opt_filter);

                    FILTERS::EXPERIMENTAL::InstAmplNormalizatorNaiveReqursive
                            <double, decltype(naive_normalizer), decltype(inst_freq_buffer),
                                    decltype(non_opt_filter)> inst_ampl_req_norm =
                            FILTERS::EXPERIMENTAL::InstAmplNormalizatorNaiveReqursive
                                    <double, decltype(naive_normalizer), decltype(inst_freq_buffer),
                                            decltype(non_opt_filter)>
                                    (naive_normalizer, inst_freq_buffer, non_opt_filter);

                    naive_normalizer.inst_freq = &inst_freq_buffer;
                    inst_ampl_req_norm.inst_freq = &inst_freq_buffer;
                    inst_ampl_req_norm.compute(*modes[iter_number], compute_buffer,
                                               &compute_buffer2);

                    for(int i = 0; i < data.size(); i++){
                        data_buffer[i] += compute_buffer[i];
                    }

                    filter.is_low_pass = true;
                    filter.compute(data_buffer, compute_buffer, compute_buffer);

                    for(int i = 0; i < data.size(); i++){
                        data_buffer[i] = data[i] - compute_buffer[i];
                        data[i] = compute_buffer[i];
                    }

                    INST_FREQ_COMPUTERS::backInstFreqNorm(data_buffer, *modes[iter_number], freq_conv);

                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);

                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);

                    inst_ampl_computer.compute(*modes[iter_number],  *inst_ampls[iter_number], nullptr);

                    iter_number++;
                }
                else{
                    INST_FREQ_COMPUTERS::backInstFreqNorm(data, *modes[iter_number], freq_conv);
                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    inst_ampl_computer.compute(*modes[iter_number],  *inst_ampls[iter_number], nullptr);

                    iter_number++;
                    break;
                }
            }
        }
    };

    struct instFreqNormSincExtractor
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

        FILTERS::InstFreqBased<double,
                FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>
                non_opt_filter = FILTERS::InstFreqBased<double,
                FILTERS::FilteringType::AverageBased,
                decltype(integrator), FILTERS::InstFreqKind::Average>(integrator);

        FILTERS::MonoFreqFilters<double, FILTERS::MonoInstFreqFilteringType::SincPaddedFIR> filter;

        template<Signal DataT>
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

            while(true){
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

                phase_computer_simple.compute(data, *phases[iter_number], nullptr);

                if((*phases[iter_number])[data.size() - 1] > 6.28){
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    double base_inst_freq = INST_FREQ_COMPUTERS::instFreqNorm(data, data_buffer, *inst_freqs[iter_number], freq_conv, freq_conv_image);

                    phase_computer_simple.compute(data_buffer, *phases[iter_number], nullptr);

                    base_inst_freq = 1.0 /
                                     (static_cast<double>(data.size()) /
                                      ((*phases[iter_number])[data.size() - 1] / 2.0 / std::numbers::pi)) / period_muller;

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

                    INST_FREQ_COMPUTERS::backInstFreqNorm(data_buffer, *modes[iter_number], freq_conv);

                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);

                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);

                    inst_ampl_computer.compute(*modes[iter_number],  *inst_ampls[iter_number], nullptr);

                    iter_number++;
                }
                else{
                    INST_FREQ_COMPUTERS::backInstFreqNorm(data, *modes[iter_number], freq_conv);

                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);
                    inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    inst_ampl_computer.compute(*modes[iter_number],  *inst_ampls[iter_number], nullptr);

                    iter_number++;
                    return;
                }
            }
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

        double period_muller = 1.1;
        double locality_coeff = 5.0;
        double max_iter_number_for_filter = 10;

        bool debug = false;

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

        FILTERS::RecursiveFilter<double, FILTERS::LocalFilteringType::SincResampled> filter;

        template<typename DataT>
        void compute(const DataT & data_in){
            filter.locality_coeff = locality_coeff;
            filter.period_muller = period_muller;
            //filter.inst_freq_computer = &inst_freq_computer;
            //filter.phase_computer = &phase_computer_simple;
            filter.max_iters = max_iter_number_for_filter;
            filter.debug = false;
            //filter.debug = debug;

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

    struct InstFreqNormSincExtractorReqV2
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

        double period_muller = 1.1;
        double locality_coeff = 5.0;
        double max_iter_number_for_filter = 10;

        bool debug = false;

        INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
        DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;

        PHASE_COMPUTERS::ExtremumsBasedUsingPCHIP
                <double, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
                phase_computer_der_atan;

        PHASE_COMPUTERS::ExtremumsBasedUsingPCHIP
                <double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
                phase_computer_simple;

        /*INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage>
                inst_freq_computer =
                INST_FREQ_COMPUTERS::ComputedOnPhase<double, decltype(integrator),
                        decltype(derivator), INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage>
                        (integrator, derivator);*/

        INST_AMPL_COMPUTERS::HilbertTransformBased
                <UTILITY_MATH::HTKind::Mull> inst_ampl_computer;

        FILTERS::SincResLocalFilterWithResReqV2<double, decltype(phase_computer_simple)> filter;

        template<typename DataT>
        void compute(const DataT & data_in){
            filter.locality_coeff = locality_coeff;
            filter.period_muller = period_muller;
            filter.phase_computer = &phase_computer_simple;
            filter.max_iters = max_iter_number_for_filter;
            filter.debug = false;
            //filter.debug = debug;

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

                if((*phases[iter_number])[data.size() - 1] >= 6.28){
                    filter.compute(data, data_buffer, &compute_buffer);

                    for(int i = 0; i < data.size(); i++){
                        (*modes[iter_number])[i] = data[i] - data_buffer[i];
                        data[i] = data_buffer[i]; //data is filtered signal
                                        //data_buffer is mode
                    }
                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);

                    //inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
                    for (size_t i = 0; i < phases[iter_number]->size(); i++){
                        (*inst_freqs[iter_number])[i] = phase_computer_simple.derive(i) / std::numbers::pi / 2.0;
                    }

                    inst_ampl_computer.compute(*modes[iter_number],  *inst_ampls[iter_number], nullptr);

                    iter_number++;
                }
                else{
                    for (auto i = 0; i < data.size(); i++){
                        (*modes[iter_number])[i] = data[i];
                    }

                    phase_computer_simple.compute(*modes[iter_number], *phases[iter_number], nullptr);
                    for (size_t i = 0; i < phases[iter_number]->size(); i++){
                        (*inst_freqs[iter_number])[i] = phase_computer_simple.derive(i) / std::numbers::pi / 2.0;
                    }
                    //inst_freq_computer.compute(*phases[iter_number], *inst_freqs[iter_number], nullptr);
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

    struct instFreqNormSincReverseSolveNonDeterminity{
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

        double period_muller = 1.1;
        double locality_coeff = 5.0;
        double max_iter_number_for_filter = 10;

        bool debug = false;

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

        FILTERS::SincResLocalFilterWithResSolveNoneDeterminityLowFreqReq
            <double, decltype(phase_computer_simple), decltype(inst_freq_computer)> filter;

        template<typename DataT>
        void compute(const DataT & data_in){
            filter.locality_coeff = locality_coeff;
            filter.period_muller = period_muller;
            filter.inst_freq_computer = &inst_freq_computer;
            filter.phase_computer = &phase_computer_simple;
            filter.max_iters = max_iter_number_for_filter;
            filter.debug = false;
            //filter.debug = debug;

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
                        /*(*modes[iter_number])[i] = data[i] - data_buffer[i];
                        data[i] = data_buffer[i]; //data is filtered signal
                                        //data_buffer is mode*/
                        (*modes[iter_number])[i] = data_buffer[i];
                        data[i] = data[i] - data_buffer[i];
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

    struct MakimaBasedModeDecomposition
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

        double max_iter_number_for_filter = 10;

        bool debug = false;

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
        
        std::vector<double> phase_shifts 
            { 0.0 * std::numbers::pi, 0.1 * std::numbers::pi, 0.2 * std::numbers::pi, 
            0.3 * std::numbers::pi, 0.4 * std::numbers::pi, 0.5 * std::numbers::pi,
            0.6 * std::numbers::pi, 0.7 * std::numbers::pi, 0.8 * std::numbers::pi, 
            0.9 * std::numbers::pi};

        FILTERS::RecursiveFilter<double, FILTERS::LocalFilteringType::MakimaInterpolationExtremums> filter;

        template<typename DataT>
        void compute(const DataT & data_in){
            filter.filter.phase_shifts = phase_shifts;
            //filter.inst_freq_computer = &inst_freq_computer;
            //filter.phase_computer = &phase_computer_simple;
            filter.max_iters = max_iter_number_for_filter;
            filter.debug = false;
            //filter.debug = debug;

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
