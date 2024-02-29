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
    template<typename U,
        Integrator<U> IntegratorT, Derivator<U> DerivatorT,
            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind kind_e>
    struct DerivativeBasedUsingExternalInstFreq {
        using AdditionalDataType = SignalPrototype<U>;

    private:
        using InstFreqDerivativeBasedKind = INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind;

    public:
        constexpr static InstFreqDerivativeBasedKind kind = kind_e;

        using BuffT = GenericSignal<SimpleVecWrapper<U>, true>;
        BuffT * inst_freq;

        IntegratorT integrator;
        DerivatorT derivator;


        DerivativeBasedUsingExternalInstFreq(IntegratorT integrator, DerivatorT derivator, BuffT& inst_freq) {
            this->integrator = integrator;
            this->derivator = derivator;
            this->inst_freq = &inst_freq;
        }

        DerivativeBasedUsingExternalInstFreq(BuffT& inst_freq) {
            this->inst_freq = &inst_freq;
        }

        DerivativeBasedUsingExternalInstFreq() {
            this->inst_freq = NULL;
        }

        template<Signal DataType, Signal OutType, Signal ComputeBufferType>
        void compute(const DataType& data, OutType& out, ComputeBufferType * computer_buffer) {
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
    template<typename U, Integrator<U> IntegratorT,
        Derivator<U> DerivatorT,
        InstFreqComputer<U> InstFreqComputerType,
        INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind kind_e>
    struct DerivativeAndInstFreqBased {
        using AdditionalDataType = SignalPrototype<U>;

    private:
        using InstFreqDerivativeBasedKind = INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind;

    public:
        constexpr static InstFreqDerivativeBasedKind kind = kind_e;

        using BuffT = GenericSignal<SimpleVecWrapper<U>, true>;
        BuffT inst_freq;

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

        template<Signal DataType, Signal OutType, Signal ComputerBufferType>
        void compute(const DataType& data, OutType& out, ComputerBufferType * computer_buffer) {
            using T = typename OutType::SampleType;
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
