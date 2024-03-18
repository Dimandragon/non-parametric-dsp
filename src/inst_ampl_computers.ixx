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

        constexpr static bool is_inst_ampl_computer = true;

        constexpr static bool is_used_external_inst_freq = true;

        using BuffT = GenericSignal<SimpleVecWrapper<U>, true>;
        BuffT * inst_freq;
        using BuffTDouble = GenericSignal<SimpleVecWrapper<std::pair<U,U>>, true>;
        BuffTDouble * inst_freq_double;

        IntegratorT integrator;
        DerivatorT derivator;


        DerivativeBasedUsingExternalInstFreq(IntegratorT integrator, DerivatorT derivator, BuffT& inst_freq) {
            this->integrator = integrator;
            this->derivator = derivator;
            this->inst_freq = &inst_freq;
        }

        DerivativeBasedUsingExternalInstFreq(IntegratorT integrator, DerivatorT derivator, BuffTDouble& inst_freq) {
            this->integrator = integrator;
            this->derivator = derivator;
            this->inst_freq_double = &inst_freq;
        }

        DerivativeBasedUsingExternalInstFreq(BuffT& inst_freq) {
            this->inst_freq = &inst_freq;
        }

        DerivativeBasedUsingExternalInstFreq(BuffTDouble& inst_freq) {
            this->inst_freq_double = &inst_freq;
        }

        DerivativeBasedUsingExternalInstFreq() {
            this->inst_freq = NULL;
        }

        template<Signal DataType, Signal OutType, Signal ComputeBufferType>
        void compute(const DataType& data, OutType& out, ComputeBufferType * computer_buffer) {
            for (int i = 0; i < out.size(); i++) {
                out[i] = 0;
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
                    out[i] = computer_buffer->interpolate(i + 0.5 / (*inst_freq)[i], SignalKind::Monotone) -
                             computer_buffer->interpolate(i - 0.5 / (*inst_freq)[i], SignalKind::Monotone) - 
                             std::abs(data.interpolate(i + 0.5 / (*inst_freq)[i], SignalKind::Monotone) - 
                                data.interpolate(i - 0.5 / (*inst_freq)[i], SignalKind::Monotone)) / 4;
                }
            } else if constexpr (kind == InstFreqDerivativeBasedKind::DeriveDouble) {
                integrator.compute(out, *computer_buffer, nullptr);
                for (int i = 0; i < out.size(); i++) {
                    out[i] = computer_buffer->interpolate(i + 0.5 / (*inst_freq_double)[i].second, SignalKind::Monotone) -
                             computer_buffer->interpolate(i - 0.5 / (*inst_freq_double)[i].first, SignalKind::Monotone) - 
                             std::abs(data.interpolate(i + 0.5 / (*inst_freq_double)[i].second, SignalKind::Monotone) - 
                                data.interpolate(i - 0.5 / (*inst_freq_double)[i].first, SignalKind::Monotone)) / 4;
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

        constexpr static bool is_inst_ampl_computer = true;
        constexpr static bool is_used_external_inst_freq = false;

        using BuffT = GenericSignal<SimpleVecWrapper<U>, true>;
        BuffT inst_freq;
        using BuffTDouble = GenericSignal<SimpleVecWrapper<std::pair<U,U>>, true>;
        BuffTDouble inst_freq_double;

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
                if constexpr (kind == InstFreqDerivativeBasedKind::DeriveDouble){
                    if (inst_freq_double.size() != data.size()) {
                        inst_freq_double.base->vec->clear();
                        for (int i = 0; i < data.size(); i++) {
                            inst_freq_double.base->vec->push_back({0.0, 0.0});
                        }
                    }
                    inst_freq_computer->compute(data, inst_freq_double, nullptr);
                }
                else{
                    if (inst_freq.size() != data.size()) {
                        static_cast<SimpleVecWrapper<T> *>(inst_freq.base)->vec->clear();
                        for (int i = 0; i < data.size(); i++) {
                            static_cast<SimpleVecWrapper<T> *>(inst_freq.base)->vec->push_back(0.);
                        }
                    }
                    inst_freq_computer->compute(data, inst_freq, nullptr);
                }
            } else {
                if constexpr (kind == InstFreqDerivativeBasedKind::DeriveDouble){
                    if (inst_freq_double.size() != data.size()) {
                        inst_freq_double.base->vec->clear();
                        for (int i = 0; i < data.size(); i++) {
                            inst_freq_double.base->vec->push_back({0.0, 0.0});
                        }
                    }
                    inst_freq_computer->compute(data, inst_freq_double, computer_buffer);
                }
                else{
                    if (inst_freq.size() != data.size()) {
                        static_cast<SimpleVecWrapper<T> *>(inst_freq.base)->vec->clear();
                        for (int i = 0; i < data.size(); i++) {
                            static_cast<SimpleVecWrapper<T> *>(inst_freq.base)->vec->push_back(0.);
                        }
                    }
                    inst_freq_computer->compute(data, inst_freq, computer_buffer);
                }
            }

            derivator.compute(data, out, nullptr);
            for (int i = 0; i < out.size(); i++) {
                out[i] = std::abs(out[i]);
            }

            if constexpr (kind == InstFreqDerivativeBasedKind::Momental
                          || kind == InstFreqDerivativeBasedKind::TimeAverage
                          || kind == InstFreqDerivativeBasedKind::DeriveAverage) {
                integrator.compute(out, *computer_buffer, nullptr);
                for (int i = 0; i < out.size(); i++) {
                    out[i] = (computer_buffer->interpolate(i + 0.5 / (inst_freq)[i], SignalKind::Monotone) -
                             computer_buffer->interpolate(i - 0.5 / (inst_freq)[i], SignalKind::Monotone) - 
                             std::abs(data.interpolate(i + 0.5 / (inst_freq)[i], SignalKind::Monotone) - 
                                data.interpolate(i - 0.5 / (inst_freq)[i], SignalKind::Monotone))) / 4;
                }
            } else if constexpr (kind == InstFreqDerivativeBasedKind::DeriveDouble) {
                integrator.compute(out, *computer_buffer, nullptr);
                for (int i = 0; i < out.size(); i++) {
                    out[i] = (computer_buffer->interpolate(i + 0.5 / (inst_freq_double)[i].second, SignalKind::Monotone) -
                             computer_buffer->interpolate(i - 0.5 / (inst_freq_double)[i].first, SignalKind::Monotone) - 
                             std::abs(data.interpolate(i + 0.5 / (inst_freq_double)[i].second, SignalKind::Monotone) - 
                                data.interpolate(i - 0.5 / (inst_freq_double)[i].first, SignalKind::Monotone))) / 4;
                }
            }
        }
    };

    export
    template<typename U, Integrator<U> Integrator, Derivator<U> Derivator, InstAmplComputer<U> InstAmplComputerT>
    struct InstAmplNormalizator{
        InstAmplComputerT * inst_ampl_computer;
        Derivator * derivator;
        Integrator * integrator;

        double min_ampl = 1.0;

        GenericSignal<SimpleVecWrapper<U>, true> buffer2;

        InstAmplNormalizator(Integrator & integrator, Derivator & derivator, InstAmplComputerT & inst_ampl_computer){
            this->inst_ampl_computer = &inst_ampl_computer;
            this->derivator = &derivator;
            this->integrator = &integrator;
        }

        template<Signal DataT, Signal OutT, Signal ComputerBufferT>
        void compute(DataT & data, OutT & out, ComputerBufferT & compute_buffer){
            auto size = data.size();
            if (buffer2.size() != size){
                buffer2.base->vec->clear();
                for(int i = 0; i < size; i++){
                    buffer2.base->vec->push_back(0);
                }
            }
            
            inst_ampl_computer->compute(data, out, &compute_buffer);

            derivator->compute(data, compute_buffer, nullptr);

            auto avg = 0.0;
            auto b = data[0];

            for (int i = 0; i < size; i++){
                avg += out[i] / size;
            }

            for (int i = 0; i < size; i++){
                if (out[i] > min_ampl){
                    compute_buffer[i] /= (out[i] / avg);
                }
            }
            integrator->compute(compute_buffer, out, nullptr);

            for (int i = 0; i < size; i++){
                out[i] += b;
            }
        }
    };
}
