module;

export module npdsp_concepts;

import <concepts>;
import <tuple>;
import <type_traits>;
import <cstddef>;
import <optional>;
import <string>;
import <utility>;
import <vector>;
import <matplot/matplot.h>;

namespace NP_DSP{
    namespace GENERAL
    {
        export
        struct Nil{};

        export
        template<typename T>
        struct Tag {
            using type = T;
        };

        template <typename T>
        constexpr bool is_any_type = true;

        export
        template <typename T>
        constexpr bool is_complex = requires(T data)
        {
            data.real();
            data.imag();
        };
    }

    namespace ONE_D{
        export enum class PlottingKind {Simple, Interpolate};

        export
        template <typename T>
        class SignalBase{
        public:
            using SampleType = T;
            using IdxType = std::size_t;
            constexpr static bool is_signal_base = true;

            virtual T & operator[](size_t idx) = 0;
            virtual T operator[](size_t idx) const = 0;
            virtual size_t size() const = 0;
            SignalBase(){}
            virtual ~SignalBase(){}
            void operator=(const SignalBase &) = delete;
            SignalBase(const SignalBase & other){std::unreachable();}
        };

        export enum class SignalKind {Monotone, Stohastic, Harmonic, Smooth};

        export
        template <typename T>
        class Signal {
        public:
            constexpr static bool is_signal = true;
            using IdxType = size_t;
            using SampleType = T;

            SignalBase<T> * base = nullptr;
            bool has_ovnership = false;

            /*inline*/
             T & operator[](size_t idx) {
                return (*base)[idx];
            }

            /*inline*/
            T operator[](size_t idx) const {
                return (*static_cast<const SignalBase<T> *>(base))[idx];
            }

            /*inline*/
            size_t size() const {
                return base->size();
            }

            //template<std::convertible_to<SignalBase<T>> BaseT>
            Signal() {
            }

            template <typename BaseT>
            Signal(GENERAL::Tag<BaseT>){
                base = new BaseT;
                has_ovnership = true;
            }

            Signal(SignalBase<T> * base) { this->base = base;  }

            virtual ~Signal() { if (has_ovnership) delete base; }
            void operator=(const Signal &) = delete;

            virtual T interpolate(double idx, SignalKind kind) const = 0;

            virtual void show(PlottingKind kind){
                if constexpr (GENERAL::is_complex<T>) {
                    if (kind == PlottingKind::Interpolate){
                        std::vector<double> plotting_data = {};
                        int i = -size();
                        while(i < static_cast<int>(size() * 2)){
                            ++i;
                            plotting_data.push_back(interpolate(static_cast<double>(i), SignalKind::Stohastic).real());
                        }
                        matplot::plot(plotting_data);
                        matplot::hold(matplot::on);
                        plotting_data.clear();
                        while(i < static_cast<int>(size() * 2)){
                            ++i;
                            plotting_data.push_back(interpolate(static_cast<double>(i), SignalKind::Stohastic).imag());
                        }
                        matplot::plot(plotting_data);
                        matplot::hold(matplot::off);
                        matplot::show();
                    }
                    else if (kind == PlottingKind::Simple){
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
                    }
                }
                else{
                    if (kind == PlottingKind::Interpolate){
                        std::vector<SampleType> plotting_data = {};
                        int i = -size();
                        while(i < static_cast<int>(size() * 2)){
                            ++i;
                            plotting_data.push_back(interpolate(static_cast<double>(i), SignalKind::Stohastic));
                        }
                        matplot::plot(plotting_data);
                        matplot::show();
                    }
                    else if (kind == PlottingKind::Simple){
                        std::vector<SampleType> plotting_data = {};
                        for (auto i = 0; i < base->size(); ++i){
                            auto sample = (*base)[i];
                            plotting_data.push_back(sample);
                        }
                        matplot::plot(plotting_data);
                        matplot::show();
                    }
                }
            }

            void show() {
                show(PlottingKind::Simple);
            }

            virtual void show(PlottingKind kind, const std::string & filename, const std::string & format) const {
                if constexpr (GENERAL::is_complex<T>) {
                    if (kind == PlottingKind::Interpolate){
                        std::vector<double> plotting_data = {};
                        int i = -size();
                        while(i < static_cast<int>(size() * 2)){
                            ++i;
                            plotting_data.push_back(interpolate(static_cast<double>(i), SignalKind::Stohastic).real());
                        }
                        matplot::plot(plotting_data);
                        matplot::hold(matplot::on);
                        plotting_data.clear();
                        while(i < static_cast<int>(size() * 2)){
                            ++i;
                            plotting_data.push_back(interpolate(static_cast<double>(i), SignalKind::Stohastic).imag());
                        }
                        matplot::plot(plotting_data);
                        matplot::hold(matplot::off);
                        matplot::show();
                        matplot::save(filename, format);
                    }
                    else if (kind == PlottingKind::Simple){
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
                        matplot::save(filename, format);
                    }
                }
                else{
                    if (kind == PlottingKind::Interpolate){
                        std::vector<SampleType> plotting_data = {};
                        int i = -size();
                        while(i < static_cast<int>(size() * 2)){
                            ++i;
                            plotting_data.push_back(interpolate(static_cast<double>(i), SignalKind::Stohastic));
                        }
                        matplot::plot(plotting_data);
                        matplot::show();
                        matplot::save(filename, format);
                    }
                    else if (kind == PlottingKind::Simple){
                        std::vector<SampleType> plotting_data = {};
                        for (auto i = 0; i < base->size(); ++i){
                            auto sample = (*base)[i];
                            plotting_data.push_back(sample);
                        }
                        matplot::plot(plotting_data);
                        matplot::show();
                        matplot::save(filename, format);
                    }
                }
            }

            virtual void show(PlottingKind kind, const std::string & filename) const {
                if constexpr (GENERAL::is_complex<T>) {
                    if (kind == PlottingKind::Interpolate){
                        std::vector<double> plotting_data = {};
                        int i = -size();
                        while(i < static_cast<int>(size() * 2)){
                            ++i;
                            plotting_data.push_back(interpolate(static_cast<double>(i), SignalKind::Stohastic).real());
                        }
                        matplot::plot(plotting_data);
                        matplot::hold(matplot::on);
                        plotting_data.clear();
                        while(i < static_cast<int>(size() * 2)){
                            ++i;
                            plotting_data.push_back(interpolate(static_cast<double>(i), SignalKind::Stohastic).imag());
                        }
                        matplot::plot(plotting_data);
                        matplot::hold(matplot::off);
                        matplot::show();
                        matplot::save(filename);
                    }
                    else if (kind == PlottingKind::Simple){
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
                        matplot::save(filename);
                    }
                }
                else{
                    if (kind == PlottingKind::Interpolate){
                        std::vector<SampleType> plotting_data = {};
                        int i = -size();
                        while(i < static_cast<int>(size() * 2)){
                            ++i;
                            plotting_data.push_back(interpolate(static_cast<double>(i), SignalKind::Stohastic));
                        }
                        matplot::plot(plotting_data);
                        matplot::show();
                        matplot::save(filename);
                    }
                    else if (kind == PlottingKind::Simple){
                        std::vector<SampleType> plotting_data = {};
                        for (auto i = 0; i < base->size(); ++i){
                            auto sample = (*base)[i];
                            plotting_data.push_back(sample);
                        }
                        matplot::plot(plotting_data);
                        matplot::show();
                        matplot::save(filename);
                    }
                }
            }

        };

        export
        template <typename T, typename SampleType>
        constexpr bool is_derivator = requires(T derivator, const Signal<SampleType> & data, Signal<SampleType> & out,
                Signal<SampleType> * additional_data)
        {
            requires T::is_derivator == true;
            typename T::DataType;
            typename T::DerivativeType;
            typename T::AdditionalDataType;

            derivator.compute(data, out, additional_data);
        };

        export
        template <typename T, typename SampleType>
        concept Derivator = is_derivator<T, SampleType>;

        export
        template <typename T, typename SampleType>
        constexpr bool is_integrator = requires(T integrator, const Signal<SampleType> & data, Signal<SampleType> & out,
                Signal<SampleType> * additional_data)
        {
            requires T::is_integrator == true;
            typename T::DataType;
            typename T::IntegralType;
            typename T::AdditionalDataType;

            integrator.compute(data, out, additional_data);
        };

        export
        template <typename T, typename SampleType>
        concept Integrator = is_integrator<T, SampleType>;

        export
        template <typename T, typename SampleType>
        constexpr bool is_modes_extracor = requires (T modes_extractor, const Signal<SampleType> & data, Signal<SampleType> & out,
                Signal<SampleType> * additional_data)
        {
            typename T::DataType;
            typename T::ModeType;
            typename T::AdditionalDataType;


            modes_extractor.compute(data, out, additional_data);
        };

        export
        template<typename T, typename SampleType>
        concept ModeExtractor = is_modes_extracor<T, SampleType>;

        export 
        template <typename T, typename SampleType>
        constexpr bool is_phase_computer
                = requires (T phase_computer, const Signal<SampleType> & data,
                    Signal<SampleType> & inst_freq, Signal<SampleType> * additional_data)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::AdditionalDataType;
            requires T::is_phase_computer == true;
            
            phase_computer.compute(data, inst_freq, additional_data);
        };

        export
        template<typename T, typename SampleType>
        concept PhaseComputer = is_phase_computer<T, SampleType>;

        export
        template <typename T, typename SampleType>
        constexpr bool is_mode_graber = requires(T mode_graber, Signal<SampleType> & data, const Signal<SampleType> & mode,
                Signal<SampleType> * additional_data)
        {
            typename T::DataType;
            typename T::ModeType;
            typename T::AdditionalDataType;
            requires T::is_mode_graber == true;

            mode_graber.compute(data, mode, additional_data);
        };

        export
        template<typename T, typename SampleType>
        concept ModeGraber = is_mode_graber<T, SampleType>;

        export
        template<typename T, typename SampleType>
        constexpr bool is_inst_freq_computer =
                requires (T inst_freq_computer, const Signal<SampleType> & data, Signal<SampleType> & inst_freq,
                Signal<SampleType> * additional_data)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::AdditionalDataType;

            requires T::is_inst_freq_computer == true;

            { inst_freq_computer.is_phase_based() } -> std::convertible_to<bool>;
            inst_freq_computer.compute(data, inst_freq, additional_data);
        };

        export
        template<typename T, typename SampleType>
        concept InstFreqComputer = is_inst_freq_computer<T, SampleType>;

        export
        template<typename T, typename SampleType>
        constexpr bool is_inst_ampl_computer =
                requires (T inst_ampl_computer, const Signal<SampleType> & data, Signal<SampleType> & inst_ampl,
                Signal<SampleType> * additional_data)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::AdditionalDataType;

            requires T::is_inst_freq_computer == true;

            inst_ampl_computer.compute(data, inst_ampl, additional_data);
        };

        export
        template<typename T, typename SampleType>
        concept InstAmplComputer = is_inst_ampl_computer<T, SampleType>;

        export
        template<typename T, typename SampleType>
        constexpr bool is_filter = requires (T filter, const Signal<SampleType> & data, Signal<SampleType> & out,
                Signal<SampleType> * additional_data)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::AdditionalDataType;

            requires T::is_filter == true;

            filter.compute(data, out, additional_data);
        };

        export
        template<typename T, typename SampleType>
        concept Filter = is_filter<T, SampleType>;


        export
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