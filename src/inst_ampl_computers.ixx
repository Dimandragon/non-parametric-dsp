module;

export module inst_ampl_computers;

import npdsp_concepts;
import utility_math;
import inst_freq_computers;
import <cmath>;
import <utility>;
import signals;
import <vector>;
import <type_traits>;
import <concepts>;

namespace NP_DSP::ONE_D::INST_AMPL_COMPUTERS {
    export
    template<typename T,
        Integrator<T> IntegratorT, Derivator <T> DerivatorT,
            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind kind_e>
    struct DerivativeBasedUsingExternalInstFreq {
        using DataType = Signal<T>;
        using OutType = Signal<T>;
        using InstFreqT = Signal<T>;
        using AdditionalDataType = Signal<T>;

    private:
        using InstFreqDerivativeBasedKind = INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind;

    public:
        constexpr static InstFreqDerivativeBasedKind kind = kind_e;

        InstFreqT* inst_freq;
        IntegratorT integrator;
        DerivatorT derivator;

        DerivativeBasedUsingExternalInstFreq(IntegratorT integrator, DerivatorT derivator, InstFreqT& inst_freq) {
            this->integrator = integrator;
            this->derivator = derivator;
            this->inst_freq = &inst_freq;
        }

        DerivativeBasedUsingExternalInstFreq(InstFreqT& inst_freq) {
            this->inst_freq = &inst_freq;
        }

        DerivativeBasedUsingExternalInstFreq() {
            this->inst_freq = NULL;
        }

        void compute(const DataType& data, OutType& out, AdditionalDataType * computer_buffer) {
            derivator.compute(data, out, nullptr);
            for (int i = 0; i < out.size(); i++) {
                out[i] = std::abs(out[i]) / 4;
            }

            if constexpr (kind == InstFreqDerivativeBasedKind::Momental
                          || kind == InstFreqDerivativeBasedKind::TimeAverage
                          || kind == InstFreqDerivativeBasedKind::DeriveAverage) {
                integrator.compute(out, *computer_buffer, nullptr);
                for (int i = 0; i < out.size(); i++) {
                    out[i] = computer_buffer->interpolate(i + 0.5 / (*inst_freq)[i]) -
                             computer_buffer->interpolate(i - 0.5 / (*inst_freq)[i]);
                }
            } else if constexpr (kind == InstFreqDerivativeBasedKind::DeriveDouble) {
                integrator.compute(out, *computer_buffer, nullptr);
                for (int i = 0; i < out.size(); i++) {
                    out[i] = computer_buffer->interpolate(i + 0.5 / (*inst_freq)[i].second) -
                             computer_buffer->interpolate(i - 0.5 / (*inst_freq)[i].first);
                }
            }
        }
    };

    export
    template<typename T, Integrator<T> IntegratorT,
        Derivator<T> DerivatorT,
        InstFreqComputer<T> InstFreqComputerType,
        INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind kind_e>
    struct DerivativeAndInstFreqBased {
        using DataType = Signal<T>;
        using OutType = Signal<T>;
        using AdditionalDataType = Signal<T>;

    private:
        using InstFreqDerivativeBasedKind = INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind;

    public:
        constexpr static InstFreqDerivativeBasedKind kind = kind_e;

        GenericSignal <T, true> inst_freq = GenericSignal <T, true>(GENERAL::Tag<SimpleVecWrapper<T>>{});

        IntegratorT integrator;
        DerivatorT derivator;
        InstFreqComputerType* inst_freq_computer;

        DerivativeAndInstFreqBased(IntegratorT integrator, DerivatorT derivator,
                                   InstFreqComputerType& inst_freq_computer) {
            this->integrator = integrator;
            this->derivator = derivator;
            this->inst_freq_computer = &inst_freq_computer;
        }

        DerivativeAndInstFreqBased(InstFreqComputerType& inst_freq_computer) {
            this->inst_freq_computer = &inst_freq_computer;
        }

        DerivativeAndInstFreqBased() {
            this->inst_freq_computer = NULL;
            //todo
        }

        //computer buffer and data must be monotone
        void compute(const DataType& data, OutType& out, AdditionalDataType * computer_buffer) {
            if constexpr (std::convertible_to<typename InstFreqComputerType::AdditionalDataType, GENERAL::Nil>) {
                if (inst_freq.size() != data.size()) {
                    static_cast<SimpleVecWrapper<T> *>(inst_freq.base)->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        static_cast<SimpleVecWrapper<T> *>(inst_freq.base)->vec->push_back(0.);
                    }
                }
                inst_freq_computer->compute(data, inst_freq, nullptr);
            } else {
                if (inst_freq.size() != data.size()) {
                    static_cast<SimpleVecWrapper<T> *>(inst_freq.base)->vec->clear();
                    for (int i = 0; i < data.size(); i++) {
                        static_cast<SimpleVecWrapper<T> *>(inst_freq.base)->vec->push_back(0.);
                    }
                }
                inst_freq_computer->compute(data, inst_freq, computer_buffer);
            }

            derivator.compute(data, out, nullptr);
            for (int i = 0; i < out.size(); i++) {
                out[i] = std::abs(out[i]) / 4;
            }

            if constexpr (kind == InstFreqDerivativeBasedKind::Momental
                          || kind == InstFreqDerivativeBasedKind::TimeAverage
                          || kind == InstFreqDerivativeBasedKind::DeriveAverage) {
                integrator.compute(out, *computer_buffer, nullptr);
                for (int i = 0; i < out.size(); i++) {
                    out[i] = computer_buffer->interpolate(i + 0.5 / inst_freq[i], SignalKind::Monotone) -
                             computer_buffer->interpolate(i - 0.5 / inst_freq[i], SignalKind::Monotone);
                }
            } else if constexpr (kind == InstFreqDerivativeBasedKind::DeriveDouble) {
                integrator.compute(out, computer_buffer, nullptr);
                for (int i = 0; i < out.size(); i++) {
                    out[i] = computer_buffer->interpolate(i + 0.5 / inst_freq[i].second, SignalKind::Monotone) -
                             computer_buffer->interpolate(i - 0.5 / inst_freq[i].first, SignalKind::Monotone);
                }
            }
        }
    };
}
