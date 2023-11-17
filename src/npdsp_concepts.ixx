module;

import <concepts>;
import <tuple>;
import <type_traits>;
import <cstddef>;

export module npdsp_concepts;

namespace NP_DSP{
    namespace GENERAL
    {
        export
        template <typename T>
        constexpr bool is_signal =
                requires (T signal, typename T::IdxType idx, std::size_t dim_number)
        {
            requires T::is_signal == true;
            requires T::dims_count == std::tuple_size_v<typename T::IdxType>;
            typename T::IdxType;
            typename T::DataType;

            { T::is_writable } -> std::same_as<bool>;

            { signal.getRefByIdx(idx) } -> std::convertible_to<typename T::DataType &>;
            { signal.getByIdx(idx) } -> std::convertible_to<typename T::DataType>;
            { signal.getDimSize(dim_number) } -> std::convertible_to<size_t>;
        };

        export
        template <typename T>
        concept Signal = is_signal<T>;
    }

    namespace ONE_D{
        export
        template <typename T>
        constexpr bool is_signal
                = requires (T signal, typename T::IdxType idx, std::size_t dim_number)
        {
            requires T::is_signal == true;
            requires std::convertible_to<typename T::IdxType, size_t>;

            typename T::IdxType;
            typename T::DataType;

            { signal.getRefByIdx(idx) } -> std::convertible_to<typename T::DataType &>;
            { signal.getByIdx(idx) } -> std::convertible_to<typename T::DataType>;
            { signal.getSize() } -> std::convertible_to<size_t>;
        };


        export
        template <typename T>
        concept Signal = is_signal<T>;
    }

    namespace GENERAL{
        export
        template <typename T>
        constexpr bool is_rotator = requires (T rotator, typename T::DataType data, typename T::OutType & out, typename T::RotorType rotor)
        {
            typename T::DataType;
            typename T::OutType;
            typename T::RotorType;

            requires  (std::tuple_size<typename T::RotorType>::value == (std::tuple_size<typename T::DataType::IdxType>::value - 1));
            requires (T::is_rotator == true);
            requires is_signal<typename T::DataType>;
            requires ! ONE_D::is_signal<typename T::DataType>;
            requires ONE_D::is_signal<typename T::OutType>;

            rotator.compute(data, out, rotor);
        };

        export
        template <typename T>
        concept Rotator = is_rotator<T>;

        export
        template <typename T>
        constexpr bool is_mode_extracor = requires (T mode_extractor, typename T::DataType data, typename T::ModeType & out)
        {
            typename T::DataType;
            typename T::ModeType;

            requires (T::is_mode_extracor == true);
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::ModeType>;

            mode_extractor.compute(data, out);
        };

        export
        template<typename T>
        concept ModeExtractor = is_mode_extracor<T>;

        export
        template <typename T>
        constexpr bool is_mode_graber = requires(T mode_graber, typename T::DataType & data, typename T::ModeType mode)
        {
            typename T::DataType;
            typename T::ModeType;

            requires T::is_mode_graber == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::ModeType>;

            mode_graber.compute(data, mode);
        };

        export
        template<typename T>
        concept ModeGraber = is_mode_graber<T>;

        export
        template<typename T>
        constexpr bool is_inst_freq_computer =
                requires (T inst_freq_computer, typename T::DataType data, typename T::OutType & inst_freq)
        {
            typename T::DataType;
            typename T::OutType;

            requires T::is_inst_freq_computer == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::InstFreqType>;

            inst_freq_computer.compute(data, inst_freq);
        };

        export
        template<typename T>
        constexpr bool is_filter = requires (T filter, typename T::DataType data, typename T::OutType & out)
        {
            typename T::DataType;
            typename T::OutType;

            requires T::is_filter == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::OutType>;

            filter.compute(data, out);
        };

        export
        template<typename T>
        concept Filter = is_filter<T>;

        export
        template<typename T>
        constexpr bool is_ortogonal_component_solver =
                requires(T solver, typename T::DataType data, typename T::InstFreqType & inst_freq,
                         typename T::InstAmplType & inst_ampl, typename T::DataType & ort_component, bool is_inst_freq_ready)
        {
            typename T::DataType;
            typename T::InstFreqType;
            typename T::InstAmplType;

            requires T::is_ortogonal_component_solver == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::InstFreqType>;
            requires is_signal<typename T::InstAmplType>;

            solver.compute(data, inst_freq, inst_ampl, ort_component, is_inst_freq_ready);
        };

        export
        template<typename T>
        concept OrtogonalComponentSolver = is_ortogonal_component_solver<T>;
    }

    namespace ONE_D
    {
        export
        template <typename T>
        constexpr bool is_derivator = requires(T derivator, typename T::DataType data, typename T::DerivativeType & out)
        {
            requires T::is_derivator == true;
            typename T::DataType;
            typename T::DerivativeType;

            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::DerivativeType>;

            derivator.compute(data, out);
        };

        export
        template <typename T>
        concept Derivator = is_derivator<T>;

        export
        template <typename T>
        constexpr bool is_integrator = requires(T integrator, typename T::DataType data, typename T::IntegralType & out)
        {
            requires T::is_integrator == true;
            typename T::DataType;
            typename T::IntegralType;

            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::IntegralType>;
            requires (std::tuple_size_v<typename T::IntegralType::IdxType> == std::tuple_size_v<typename T::DataType::IdxType>);

            integrator.compute(data, out);
        };

        export
        template <typename T>
        concept Integrator = is_integrator<T>;

        export
        template <typename T>
        constexpr bool is_mode_extracor = requires (T mode_extractor, typename T::DataType data, typename T::ModeType & out)
        {
            typename T::DataType;
            typename T::ModeType;

            requires (T::is_mode_extracor == true);
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::ModeType>;

            mode_extractor.compute(data, out);
        };

        export
        template<typename T>
        concept ModeExtractor = is_mode_extracor<T>;

        export
        template <typename T>
        constexpr bool is_mode_graber = requires(T mode_graber, typename T::DataType & data, typename T::ModeType mode)
        {
            typename T::DataType;
            typename T::ModeType;

            requires T::is_mode_graber == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::ModeType>;

            mode_graber.compute(data, mode);
        };

        export
        template<typename T>
        concept ModeGraber = is_mode_graber<T>;

        export
        template<typename T>
        constexpr bool is_inst_freq_computer =
                requires (T inst_freq_computer, typename T::DataType data, typename T::OutType & inst_freq)
        {
            typename T::DataType;
            typename T::OutType;

            requires T::is_inst_freq_computer == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::InstFreqType>;

            inst_freq_computer.compute(data, inst_freq);
        };

        export
        template<typename T>
        concept InstFreqComputer = is_inst_freq_computer<T>;

        export
        template<typename T>
        constexpr bool is_filter = requires (T filter, typename T::DataType data, typename T::OutType & out)
        {
            typename T::DataType;
            typename T::OutType;

            requires T::is_filter == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::OutType>;

            filter.compute(data, out);
        };

        export
        template<typename T>
        concept Filter = is_filter<T>;

        export
        template<typename T>
        constexpr bool is_ortogonal_component_solver =
                requires(T solver, typename T::DataType data, typename T::InstFreqType & inst_freq,
                        typename T::InstAmplType & inst_ampl, typename T::DataType & ort_component, bool is_inst_freq_ready)
        {
            typename T::DataType;
            typename T::InstFreqType;
            typename T::InstAmplType;

            requires T::is_ortogonal_component_solver == true;
            requires is_signal<typename T::DataType>;
            requires is_signal<typename T::InstFreqType>;
            requires is_signal<typename T::InstAmplType>;

            solver.compute(data, inst_freq, inst_ampl, ort_component, is_inst_freq_ready);
        };

        export
        template<typename T>
        concept OrtogonalComponentSolver = is_ortogonal_component_solver<T>;
    }
}


//static_assert(is_signal<vec_wrapper<int>, int, int>);