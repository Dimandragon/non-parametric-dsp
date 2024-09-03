#pragma once

#include <concepts>
#include <type_traits>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace NP_DSP{
    namespace GENERAL
    {
        
        struct Nil{};

        
        template<typename T>
        struct Tag {
            using type = T;
        };

        template <typename T>
        constexpr bool is_any_type = true;

        
        template <typename T>
        constexpr bool is_complex = requires(T data)
        {
            data.real();
            data.imag();
        };
    }

    namespace ONE_D{
        
        template <typename T>
        constexpr bool is_signal_base_first
                = requires (T signal, T::IdxType idx)
        {
            requires T::is_signal_base == true;
            //requires T::is_writable == true || T::is_writable == false;

            typename T::IdxType;
            typename T::SampleType ;

            //{ delete new T };
            //{ signal[idx] } -> std::convertible_to<typename T::SampleType &>;


            { signal.size() } -> std::convertible_to<size_t>;
        };

        
        template <typename T>
        constexpr bool is_signal_base_second
                = requires (T signal, T::IdxType idx)
        {
            requires T::is_signal_base == true;
            requires (T::is_writable == true || T::is_writable == false);

            typename T::IdxType;
            typename T::SampleType;

            { delete new T };
            { signal[idx] } -> std::convertible_to<typename T::SampleType>;

            { signal.size() } -> std::convertible_to<size_t>;
        };



        
        template <typename T>
        concept SignalBase = is_signal_base_second<T> || is_signal_base_first<T>;

         enum class SignalKind {Monotone, Stohastic, Harmonic, Smooth, Universal};
         enum class PlottingKind {Simple, Interpolate, Spectre};

        
        template <typename T>
        constexpr bool is_signal = requires(T signal, size_t idx, std::optional<typename T::IdxType> idx1,
        T::SampleType value, PlottingKind kind, const std::string & filename, const std::string & format, SignalKind s_kind) {
            requires T::is_signal == true;

            //typename T::Base;
            //requires is_signal_base_first<typename T::Base> || is_signal_base_second<typename T::Base>;

            typename T::IdxType;
            requires std::is_same_v<typename T::IdxType, typename T::Base::IdxType>;

            typename T::SampleType;
            requires std::is_same_v<typename T::SampleType, typename T::Base::SampleType>;

            //{ delete new T };
            //{ signal[idx] } -> std::convertible_to<typename T::SampleType &>;
            { signal.size() } -> std::convertible_to<size_t>;
            { signal.interpolate(idx, s_kind) } -> std::convertible_to<typename T::SampleType>;
            //{ signal.findInterpolateUnimode(value, idx, idx) } -> std::convertible_to<typename T::IdxType>;
            //{ signal.findMonotone(value, idx1, idx1)} -> std::convertible_to<typename T::IdxType>;
            { signal.show(kind, filename) };
            { signal.show(kind, filename, format) };
        };

        
        template<typename T>
        concept Signal = is_signal<T>;

        
        template <typename T>
        constexpr bool is_signal_wrapper = requires{
            requires is_signal<T> || std::is_same_v<GENERAL::Nil, T>;
        };

        
        template <typename T>
        concept SignalWrapper = is_signal_wrapper<T>;

        
        template <typename T>
        class SignalBasePrototype{
        public:
            using SampleType = T;
            using IdxType = size_t;
            constexpr static bool is_signal_base = true;

            virtual T & operator[](size_t idx) = 0;
            virtual T operator[](size_t idx) const = 0;
            virtual size_t size() const = 0;
            SignalBasePrototype(){}
            virtual ~SignalBasePrototype(){}
            void operator=(const SignalBasePrototype &) = delete;
            //SignalBasePrototype(const SignalBasePrototype & other){/*std::unreachable();*/}
        };

        
        template <typename T>
        class SignalPrototype {
        public:
            constexpr static bool is_signal = true;
            using IdxType = size_t;
            using SampleType = T;
            using Base = SignalBasePrototype<T>;

            SignalBasePrototype<T> * base = nullptr;
            bool has_ovnership = false;

            /*inline*/
             T & operator[](size_t idx) {
                return (*base)[idx];
            }

            /*inline*/
            T operator[](size_t idx) const {
                return (*static_cast<const SignalBasePrototype<T> *>(base))[idx];
            }

            /*inline*/
            size_t size() const {
                return base->size();
            }

            //template<std::convertible_to<SignalBase<T>> BaseT>
            SignalPrototype() {
            }

            template <typename BaseT>
            SignalPrototype(GENERAL::Tag<BaseT>){
                base = new BaseT;
                base = new BaseT;
                has_ovnership = true;
            }

            SignalPrototype(SignalBasePrototype<T> * base) { this->base = base;  }

            virtual ~SignalPrototype() { if (has_ovnership) delete base; }
            void operator=(const SignalPrototype &) = delete;

            virtual T interpolate(double idx, SignalKind kind) const = 0;

            virtual void show(PlottingKind kind){
                if constexpr (GENERAL::is_complex<T>) {
                    if (kind == PlottingKind::Interpolate){
                        
                    }
                    /*else if (kind == PlottingKind::Simple){
                        std::vector<double> plotting_data = {};
                        for (auto i = 0; i < base->size(); ++i){
                            auto sample = (*base)[i];
                            plotting_data.push_back(sample.real());
                        }
                        matplot::plot(plotting_data);
                        matplot::hold(matplot::on);
                        plotting_data.clear();
                        for (auto i = 0; i < base->size(); ++i){
                            auto sample = (*base)[i];
                            plotting_data.push_back(sample.imag());
                        }
                        matplot::plot(plotting_data);
                        matplot::hold(matplot::off);
                        matplot::show();
                    }*/
                }
                else{
                    if (kind == PlottingKind::Interpolate){
                        
                    }
                    else if (kind == PlottingKind::Simple){
                        
                    }
                }
            }

            void show() {
                show(PlottingKind::Simple);
            }

            virtual void show(PlottingKind kind, const std::string & filename, const std::string & format) const {
                if constexpr (GENERAL::is_complex<T>) {
                    if (kind == PlottingKind::Interpolate){
                        
                    }
                    else if (kind == PlottingKind::Simple){
                        
                    }
                }
                else{
                    if (kind == PlottingKind::Interpolate){
                        
                    }
                    else if (kind == PlottingKind::Simple){
                        
                    }
                }
            }

            virtual void show(PlottingKind kind, const std::string & filename) const {
                if constexpr (GENERAL::is_complex<T>) {
                    if (kind == PlottingKind::Interpolate){
                        
                    }
                    else if (kind == PlottingKind::Simple){
                        
                    }
                }
                else{
                    if (kind == PlottingKind::Interpolate){
                        
                    }
                    else if (kind == PlottingKind::Simple){
                        
                    }
                }
            }
        };
        struct details {
            static_assert(is_signal<SignalPrototype<double>>);
        };

        
        template <typename T, typename SampleType>
        constexpr bool is_derivator = requires(T derivator, const SignalPrototype<SampleType> & data, SignalPrototype<SampleType> & out,
                SignalPrototype<SampleType> * additional_data)
        {
            typename T::AdditionalDataType;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            requires T::is_derivator == true;

            derivator.compute(data, out, additional_data);
        };

        
        template <typename T, typename SampleType>
        concept Derivator = is_derivator<T, SampleType>;

        
        template <typename T, typename SampleType>
        constexpr bool is_integrator = requires(T integrator, const SignalPrototype<SampleType> & data, SignalPrototype<SampleType> & out,
                SignalPrototype<SampleType> * additional_data)
        {
            requires T::is_integrator == true;
            typename T::AdditionalDataType;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            integrator.compute(data, out, additional_data);
        };

        
        template <typename T, typename SampleType>
        concept Integrator = is_integrator<T, SampleType>;

        
        template <typename T, typename SampleType>
        constexpr bool is_modes_extracor = requires (T modes_extractor, const SignalPrototype<SampleType> & data, SignalPrototype<SampleType> & out,
                SignalPrototype<SampleType> * additional_data)
        {
            requires T::is_modes_extractor == true;
            typename T::AdditionalDataType;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            modes_extractor.compute(data, out, additional_data);
        };

        
        template<typename T, typename SampleType>
        concept ModeExtractor = is_modes_extracor<T, SampleType>;

         
        template <typename T, typename SampleType>
        constexpr bool is_phase_computer
                = requires (T phase_computer, const SignalPrototype<SampleType> & data,
                    SignalPrototype<SampleType> & inst_freq, SignalPrototype<SampleType> * additional_data)
        {
            requires T::is_phase_computer == true;
            typename T::AdditionalDataType;
            requires is_signal_wrapper<typename T::AdditionalDataType>;
            
            phase_computer.compute(data, inst_freq, additional_data);
        };

        
        template<typename T, typename SampleType>
        concept PhaseComputer = is_phase_computer<T, SampleType>;

        
        template <typename T, typename SampleType>
        constexpr bool is_mode_graber = requires(T mode_graber, SignalPrototype<SampleType> & data, const SignalPrototype<SampleType> & mode,
                SignalPrototype<SampleType> * additional_data)
        {
            requires T::is_mode_graber == true;
            typename T::AdditionalDataType;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            mode_graber.compute(data, mode, additional_data);
        };

        
        template<typename T, typename SampleType>
        concept ModeGraber = is_mode_graber<T, SampleType>;

        
        template<typename T, typename SampleType>
        constexpr bool is_inst_freq_computer =
                requires (T inst_freq_computer, const SignalPrototype<SampleType> & data, SignalPrototype<SampleType> & inst_freq,
                SignalPrototype<SampleType> * additional_data)
        {
            requires T::is_inst_freq_computer == true;
            typename T::AdditionalDataType;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            { inst_freq_computer.is_phase_based() } -> std::convertible_to<bool>;
            inst_freq_computer.compute(data, inst_freq, additional_data);
        };

        
        template<typename T, typename SampleType>
        concept InstFreqComputer = is_inst_freq_computer<T, SampleType>;

        
        template<typename T, typename SampleType>
        constexpr bool is_inst_ampl_computer =
                requires (T inst_ampl_computer, const SignalPrototype<SampleType> & data, SignalPrototype<SampleType> & inst_ampl,
                SignalPrototype<SampleType> * additional_data)
        {
            requires T::is_inst_ampl_computer == true;
            typename T::AdditionalDataType;

            inst_ampl_computer.compute(data, inst_ampl, additional_data);
        };

        
        template<typename T, typename SampleType>
        concept InstAmplComputer = is_inst_ampl_computer<T, SampleType>;

        
        template<typename T, typename SampleType>
        constexpr bool is_filter = requires (T filter, const SignalPrototype<SampleType> & data, SignalPrototype<SampleType> & out,
                SignalPrototype<SampleType> * additional_data)
        {
            requires T::is_filter == true;
            typename T::AdditionalDataType;

            filter.compute(data, out, additional_data);
        };

        
        template<typename T, typename SampleType>
        concept Filter = is_filter<T, SampleType>;


        
        template<typename T, typename SampleType>
        constexpr bool is_signal_approximator =
                requires(T & approximator, T::Loss loss, T::StopPoint stop_point, T::ApproxModel model, T::IdxType idx, T::SampleType max_error){
            typename T::IdxType;
            typename T::SampleType;
            typename T::SignalType;
            typename T::Loss;
            typename T::StopPoint;
            //typename T::ApproxModel;

            requires T::is_signal_approximator == true;

            { approximator.latest_loss } -> std::convertible_to<typename T::SampleType>;
            { approximator.max_error } -> std::convertible_to<typename T::SampleType>;
            { delete new T(loss, stop_point, model, max_error) };
            { stop_point(loss, approximator) } -> std::convertible_to<bool>;
            //here loss is the opt iteration losses different
            { loss(approximator) } -> std::convertible_to<typename T::SampleType>;
            { approximator.train() };
            { approximator.compute(idx) } -> std::convertible_to<typename T::SampleType>;
        };
    }
}


//static_assert(is_signal<vec_wrapper<int>, int, int>);