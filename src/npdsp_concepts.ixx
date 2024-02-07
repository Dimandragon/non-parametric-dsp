module;

export module npdsp_concepts;

import <concepts>;
import <tuple>;
import <type_traits>;
import <cstddef>;
import <optional>;
import <string>;

namespace NP_DSP{
    namespace GENERAL
    {
        export
        struct Nil{};

        export
        template <typename T>
        constexpr bool is_signal_base =
                requires (T signal, typename T::IdxType idx, std::size_t dim_number)
        {
            requires (T::is_writable == true || T::is_writable == false);
            requires T::is_signal == true;
            requires T::dims_count == std::tuple_size_v<typename T::IdxType>;
            typename T::IdxType;
            typename T::SampleType;

            { delete new T };
            { signal[idx] } -> std::convertible_to<typename T::SampleType &>;
            { signal.dimSize(dim_number) } -> std::convertible_to<size_t>;
        };

        export
        template <typename T>
        concept SignalBase = is_signal_base<T>;

        export enum class PlottingKind {Simple, Interpolate};

        export
        template <typename T>
        constexpr bool is_signal = requires (T signal, T::IdxType idx, size_t dim_number,
                std::optional<typename T::IdxType> idx1, std::optional<typename T::IdxType> idx2,
                T::SampelType  value, PlottingKind kind, const std::string & filename, const std::string & format)
        {
            requires T::is_signal == true;

            typename T::Base;
            requires is_signal_base<typename T::Base>;

            typename T::IdxType;
            requires std::is_same_v<typename T::IdxType, typename T::Base::IdxType>;

            typename T::SampleType;
            requires std::is_same_v<typename T::SampleType, typename T::Base::IdxType>;

            { delete new T };
            { signal[idx] } -> std::convertible_to<typename T::SampleType &>;
            { signal.dimSize(dim_number) } -> std::convertible_to<size_t>;
            { signal.interpolate(idx) } -> std::convertible_to<typename T::SampleType>;
            { signal.findInterpolate(value, idx1, idx2) } -> std::convertible_to<typename T::IdxType>;
            { signal.findMonotone(value, idx1, idx2)} -> std::convertible_to<typename T::IdxType>;
            signal.show(kind, filename);
            signal.show(kind, filename, format);
        };

        export
        template<typename T>
        concept Signal = is_signal<T>;

        export
        template <typename T>
        constexpr bool is_signal_wrapper = requires{
            requires is_signal<T> || std::is_same_v<GENERAL::Nil, T>;
        };

        export
        template <typename T>
        concept SignalWrapper = is_signal_wrapper<T>;
    }

    namespace ONE_D{
        export
        template <typename T>
        constexpr bool is_signal_base_first
                = requires (T signal, T::IdxType idx)
        {
            requires T::is_signal_base == true;
            requires (T::is_writable == true || T::is_writable == false);

            typename T::IdxType;
            typename T::SampleType;

            //{ delete new T };
            //{ signal[idx] } -> std::convertible_to<typename T::SampleType &>;
            
            
            { signal.size() } -> std::convertible_to<size_t>;
        };

        export
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

        

        export
        template <typename T>
        concept SignalBase = (is_signal_base_second<T> || is_signal_base_first<T>);

        export enum class PlottingKind {Simple, Interpolate};

        export
        template <typename T>
        constexpr bool is_signal = requires(T signal, T::IdxType idx, std::optional<typename T::IdxType> idx1,
        T::SampleType value, PlottingKind kind, const std::string & filename, const std::string & format) {
            requires T::is_signal == true;

            typename T::Base;
            requires (is_signal_base_first<typename T::Base> || is_signal_base_second<typename T::Base>);

            typename T::IdxType;
            requires std::is_same_v<typename T::IdxType, typename T::Base::IdxType>;

            typename T::SampleType;
            requires std::is_same_v<typename T::SampleType, typename T::Base::SampleType>;

            { delete new T };
            //{ signal[idx] } -> std::convertible_to<typename T::SampleType &>;
            { signal.size() } -> std::convertible_to<size_t>;
            { signal.interpolate(idx) } -> std::convertible_to<typename T::SampleType>;
            //{ signal.findInterpolateUnimode(value, idx, idx) } -> std::convertible_to<typename T::IdxType>;
            { signal.findMonotone(value, idx1, idx1)} -> std::convertible_to<typename T::IdxType>;
            { signal.show(kind, filename) };
            { signal.show(kind, filename, format) };
        };  

        export
        template<typename T>
        concept Signal = is_signal<T>;


        export
        template <typename T>
        constexpr bool is_signal_wrapper = requires{
            requires is_signal<T> || std::is_same_v<GENERAL::Nil, T>;
        };

        export
        template <typename T>
        concept SignalWrapper = is_signal_wrapper<T>;
    }

    namespace GENERAL{
        export
        template <typename T>
        constexpr bool is_rotator =
                requires (T rotator, const T::DataType & data, T::OutType & out, T::RotorType rotor, T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::RotorType;
            typename T::AdditionalDataType;

            requires  (std::tuple_size<typename T::RotorType>::value == (std::tuple_size<typename T::DataType::IdxType>::value - 1));
            requires (T::is_rotator == true);
            requires is_signal<typename T::DataType>;
            requires ! ONE_D::is_signal<typename T::DataType>;
            requires ONE_D::is_signal<typename T::OutType>;
            requires (T::OutType::is_writable == true);
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            rotator.compute(data, out, rotor, additional_data);
        };

        export
        template <typename T>
        concept Rotator = is_rotator<T>;

        export
        template <typename T>
        constexpr bool is_mode_extracor =
                requires (T mode_extractor, const T::DataType & data, T::ModeType & out, T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::ModeType;
            typename T::AdditionalDataType;

            requires (T::is_mode_extracor == true);
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::ModeType>;
            requires is_signal_wrapper<typename T::AdditionalDataType>;
            requires (T::ModeType::is_writable == true);

            mode_extractor.compute(data, out, additional_data);
        };

        export
        template<typename T>
        concept ModeExtractor = is_mode_extracor<T>;

        export
        template <typename T>
        constexpr bool is_mode_graber =
                requires(T mode_graber, T::DataType & data, const T::ModeType & mode, T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::ModeType;
            typename T::AdditionalDataType;

            requires T::is_mode_graber == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::ModeType>;
            requires (T::DataType::is_writable == true);
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            mode_graber.compute(data, mode, additional_data);
        };

        export
        template<typename T>
        concept ModeGraber = is_mode_graber<T>;

        export
        template<typename T>
        constexpr bool is_inst_freq_computer =
                requires (T inst_freq_computer, const T::DataType & data, T::OutType & inst_freq, T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::AdditionalDataType;

            requires T::is_inst_freq_computer == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::InstFreqType>;
            requires (T::OutType::is_writable == true);
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            inst_freq_computer.compute(data, inst_freq, additional_data);
        };

        export
        template<typename T>
        constexpr bool is_filter = requires (T filter, const T::DataType & data, T::OutType & out, T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::InstFreqType;
            typename T::AdditionalDataType;

            requires T::is_filter == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::OutType>;
            requires is_signal<typename T::InstFreqType>;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            filter.compute(data, out, additional_data);
        };

        export
        template<typename T>
        concept Filter = is_filter<T>;

        export
        template<typename T>
        constexpr bool is_ortogonal_component_solver =
                requires(T solver, const T::DataType & data, T::InstFreqType & inst_freq, T::InstAmplType & inst_ampl,
                        T::DataType & ort_component, bool is_inst_freq_ready, T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::InstFreqType;
            typename T::InstAmplType;
            typename T::AdditionalDataType;

            requires T::is_ortogonal_component_solver == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::InstFreqType>;
            requires is_signal<typename T::InstAmplType>;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            solver.compute(data, inst_freq, inst_ampl, ort_component, is_inst_freq_ready, additional_data);
        };

        export
        template<typename T>
        concept OrtogonalComponentSolver = is_ortogonal_component_solver<T>;
    }

    namespace ONE_D
    {
        export
        template <typename T>
        constexpr bool is_derivator = requires(T derivator, const T::DataType & data, T::DerivativeType & out,
                T::AdditionalDataType & additional_data)
        {
            requires T::is_derivator == true;
            typename T::DataType;
            typename T::DerivativeType;
            typename T::AdditionalDataType;

            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::DerivativeType>;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            derivator.compute(data, out, additional_data);
        };

        export
        template <typename T>
        concept Derivator = is_derivator<T>;

        export
        template <typename T>
        constexpr bool is_integrator = requires(T integrator, const T::DataType & data, T::IntegralType & out,
                T::AdditionalDataType & additional_data)
        {
            requires T::is_integrator == true;
            typename T::DataType;
            typename T::IntegralType;
            typename T::AdditionalDataType;

            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::IntegralType>;
            //requires (std::tuple_size_v<typename T::IntegralType::IdxType> == std::tuple_size_v<typename T::DataType::IdxType>);
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            integrator.compute(data, out, additional_data);
        };

        export
        template <typename T>
        concept Integrator = is_integrator<T>;

        export
        template <typename T>
        constexpr bool is_mode_extracor = requires (T mode_extractor, const T::DataType & data, T::ModeType & out,
                T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::ModeType;
            typename T::AdditionalDataType;

            requires (T::is_mode_extracor == true);
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::ModeType>;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            mode_extractor.compute(data, out, additional_data);
        };

        export
        template<typename T>
        concept ModeExtractor = is_mode_extracor<T>;

        export
        template <typename T>
        constexpr bool is_mode_graber = requires(T mode_graber, T::DataType & data, const T::ModeType & mode,
                T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::ModeType;
            typename T::AdditionalDataType;
            requires T::is_mode_graber == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::ModeType>;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            mode_graber.compute(data, mode, additional_data);
        };

        export
        template<typename T>
        concept ModeGraber = is_mode_graber<T>;

        export
        template<typename T>
        constexpr bool is_inst_freq_computer =
                requires (T inst_freq_computer, const T::DataType & data, T::OutType & inst_freq,
                T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::AdditionalDataType;

            requires T::is_inst_freq_computer == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::InstFreqType>;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            inst_freq_computer.compute(data, inst_freq, additional_data);
        };

        export
        template<typename T>
        concept InstFreqComputer = is_inst_freq_computer<T>;

        export
        template<typename T>
        constexpr bool is_filter = requires (T filter, const T::DataType & data, T::OutType & out,
                T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::AdditionalDataType;

            requires T::is_filter == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::OutType>;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            filter.compute(data, out, additional_data);
        };

        export
        template<typename T>
        concept Filter = is_filter<T>;

        export
        template<typename T>
        constexpr bool is_ortogonal_component_solver =
                requires(T solver, const T::DataType & data, T::InstFreqType & inst_freq,
                        T::InstAmplType & inst_ampl, T::DataType & ort_component, bool is_inst_freq_ready,
                        T::AdditionalDataType & additional_data)
        {
            typename T::DataType;
            typename T::InstFreqType;
            typename T::InstAmplType;
            typename T::AdditionalDataType;

            requires T::is_ortogonal_component_solver == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::InstFreqType>;
            requires is_signal<typename T::InstAmplType>;
            requires is_signal_wrapper<typename T::AdditionalDataType>;

            solver.compute(data, inst_freq, inst_ampl, ort_component, is_inst_freq_ready, additional_data);
        };

        export
        template<typename T>
        concept OrtogonalComponentSolver = is_ortogonal_component_solver<T>;

        export
        template<typename T>
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