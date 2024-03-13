module;

#include <icecream.hpp>

export module filters;

import npdsp_concepts;
import signals;
import <utility>;
import integrators;
import <string>;
import inst_freq_computers;
import utility_math;
import <concepts>;
import phase_computers;
import <vector>;
import <algorithm>;
import inst_ampl_computers;


namespace NP_DSP::ONE_D::FILTERS {
    export enum class InstFreqKind { Average, Double };

    export enum class FilteringType { DerivativeBased, ValueBased, AverageBased, Median, ValueBasedSmart, DerivativeBasedSmart };

    export
    template<typename U, FilteringType filtering_type_k,
        Integrator<U> IntegratorT, InstFreqKind inst_freq_k>
    struct NonOptPeriodBasedFilter {
        using AdditionalDataType = SignalPrototype<double>;
        using IdxType = size_t;

        IntegratorT integrator;
        constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
        constexpr static bool is_filter = true;

        double step = 0.1;

        NonOptPeriodBasedFilter(IntegratorT integrator) {
            this->integrator = integrator;
        }

        NonOptPeriodBasedFilter() {
        }

        //data and inst freq must be not monotone
        template<Signal DataType, Signal OutType, Signal InstFreqType>
        void compute(const DataType& data, OutType& out, const InstFreqType * inst_freq) {
            using T = typename OutType::SampleType;
            if constexpr (filtering_type_k == FilteringType::DerivativeBased) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    auto size_expr = [&]() {
                        return data.size();
                    };

                    auto val_expression = [&](size_t idx) {
                        //IC(inst_freq);
                        return (data.interpolate(idx + 0.5 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic) - data.interpolate(
                                    idx - 0.5 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) *
                               inst_freq->interpolate(idx, SignalKind::Stohastic);
                    };

                    ExpressionWrapper<T, size_t,
                        decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                            expr_wrapper(val_expression, size_expr);

                    const GenericSignal<decltype(expr_wrapper), false> expr_signal(expr_wrapper);

                    INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator_new;
                    //todo fix its big crucher
                    integrator_new.compute(expr_signal, out, nullptr);

                    double avg = 0.0;
                    for (int i = 0; i < data.size(); i++){
                        avg += out[i];
                    }
                    avg = avg / data.size();

                    for (int i = 0; i < data.size(); i++){
                        out[i] -= avg;
                    }
                } else if constexpr (inst_freq_kind == InstFreqKind::Double) {
                    auto inst_freq_first_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].first;
                    };
                    auto inst_freq_second_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].second;
                    };
                    auto size_expr = [&]() {
                        return data.size();
                    };
                    ExpressionWrapper<T, size_t,
                        decltype(inst_freq_first_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_first_expr_wrapper(inst_freq_first_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_first_expr_wrapper), false> inst_freq_first(inst_freq_first_expr_wrapper);

                    ExpressionWrapper<T, size_t,
                        decltype(inst_freq_second_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_second_expr_wrapper(inst_freq_second_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_second_expr_wrapper), false> inst_freq_second (inst_freq_second_expr_wrapper);

                    auto val_expression = [&](typename DataType::IdxType idx) {
                        return (data.interpolate(idx + 0.5 / inst_freq_second.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic) -
                                data.interpolate(idx - 0.5 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) /
                               (0.5 / inst_freq_second.interpolate(idx, SignalKind::Stohastic) + 0.5 / inst_freq_first.interpolate(idx, SignalKind::Stohastic));
                    };

                    ExpressionWrapper<T, size_t,
                                decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                            expr_wrapper(val_expression, size_expr);

                    GenericSignal<decltype(expr_wrapper), false> expr_signal(expr_wrapper);

                    INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator_new;

                    //todo fix its big crucher
                    integrator_new.compute(expr_signal, out, nullptr);

                    double avg = 0.0;
                    for (int i = 0; i < data.size(); i++){
                        avg += out[i];
                    }
                    avg = avg / data.size();

                    for (int i = 0; i < data.size(); i++){
                        out[i] -= avg;
                    }
                } else {
                    std::unreachable();
                }
            } else if constexpr (filtering_type_k == FilteringType::ValueBased) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    auto val_expression = [&](typename DataType::IdxType idx) {
                        return (data.interpolate(idx + 0.25 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)
                                + data.interpolate(idx - 0.25 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) / 2;
                    };
                    for (auto i = 0; i < data.size(); i++) {
                        out[i] = val_expression(i);
                    }
                } else if constexpr (inst_freq_kind == InstFreqKind::Double) {
                    auto inst_freq_first_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].first;
                    };
                    auto inst_freq_second_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].second;
                    };
                    auto size_expr = [&]() {
                        return data.size();
                    };
                    ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_first_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_first_expr_wrapper(inst_freq_first_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_first_expr_wrapper), false> inst_freq_first(
                        inst_freq_first_expr_wrapper);

                    ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_second_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_second_expr_wrapper(inst_freq_second_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_second_expr_wrapper), false> inst_freq_second(
                        inst_freq_second_expr_wrapper);

                    auto val_expression = [&](typename DataType::IdxType idx) {
                        return (data.interpolate(idx + 0.25 / inst_freq_second.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic) +
                                data.interpolate(idx - 0.25 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) / 2.0;
                        /*return UTILITY_MATH::linearInterpolate<double, double>
                        ({idx - 0.25 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), 
                            data.interpolate(idx - 0.25 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)},
                            {idx + 0.25 / inst_freq_second.interpolate(idx, SignalKind::Stohastic),
                            data.interpolate(idx + 0.25 / inst_freq_second.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)},
                            idx);*/
                    };
                    for (auto i = 0; i < data.size(); i++) {
                        out[i] = val_expression(i);
                    }
                    
                } else {
                    std::unreachable();
                }
            } else if constexpr (filtering_type_k == FilteringType::ValueBasedSmart) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    auto val_expression = [&](typename DataType::IdxType idx) {
                        for (int i = idx; i >= 0; i--){
                            //find 
                            if (i + 0.25 / inst_freq->interpolate(i, SignalKind::Stohastic) <= idx){
                                std::pair<double, double> point1 = {static_cast<double>(i), i + 0.25 / inst_freq->interpolate(i, SignalKind::Stohastic)};
                                std::pair<double, double> point2 = {static_cast<double>(i+1), i + 1. + 0.25 / inst_freq->interpolate(i+1, SignalKind::Stohastic)};
                                auto dy = (point2.second - point1.second);
                                auto dx = (point2.first - point1.first);
                                double idx_new;

                                if (dy == 0) {
                                    idx_new = (point1.first + point2.first) / 2;
                                }
                                idx_new = point1.first + dx * (idx - point1.second) / dy;
                                
                                auto dy1_new = data.interpolate(idx_new, SignalKind::Stohastic) - 
                                    data.interpolate(idx_new - 0.25 / inst_freq->interpolate(idx_new, SignalKind::Stohastic), SignalKind::Stohastic);
                                auto dx1_new = 0.25 / inst_freq->interpolate(idx_new, SignalKind::Stohastic);

                                auto dx2_new = 0.25 / inst_freq->interpolate(idx_new, SignalKind::Stohastic);
                                auto dy2_new = dy1_new;
                                return data.interpolate(idx_new, SignalKind::Stohastic) + dy2_new / dx1_new * dx2_new ;
                            }
                        }
                        
                        
                        return (data.interpolate(idx + 0.25 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)
                                + data.interpolate(idx - 0.25 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) / 2;
                    };
                    for (auto i = 0; i < data.size(); i++) {
                        out[i] = val_expression(i);
                    }
                } else if constexpr (inst_freq_kind == InstFreqKind::Double) {
                    auto inst_freq_first_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].first;
                    };
                    auto inst_freq_second_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].second;
                    };
                    auto size_expr = [&]() {
                        return data.size();
                    };
                    ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_first_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_first_expr_wrapper(inst_freq_first_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_first_expr_wrapper), false> inst_freq_first(
                        inst_freq_first_expr_wrapper);

                    ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_second_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_second_expr_wrapper(inst_freq_second_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_second_expr_wrapper), false> inst_freq_second(
                        inst_freq_second_expr_wrapper);

                    auto val_expression = [&](typename DataType::IdxType idx) {
                        for (int i = idx; i >= 0; i--){
                            //find 
                            if (i + 0.25 / inst_freq_second.interpolate(i, SignalKind::Stohastic) <= idx){
                                std::pair<double, double> point1 = {static_cast<double>(i), i + 0.25 / inst_freq_second.interpolate(i, SignalKind::Stohastic)};
                                std::pair<double, double> point2 = {static_cast<double>(i+1), i + 1. + 0.25 / inst_freq_second.interpolate(i+1, SignalKind::Stohastic)};
                                auto dy = (point2.second - point1.second);
                                auto dx = (point2.first - point1.first);
                                double idx_new;

                                if (dy == 0) {
                                    idx_new = (point1.first + point2.first) / 2;
                                }
                                idx_new = point1.first + dx * (idx - point1.second) / dy;
                                
                                auto dy1_new = data.interpolate(idx_new, SignalKind::Stohastic) - 
                                    data.interpolate(idx_new - 0.25 / inst_freq_first.interpolate(idx_new, SignalKind::Stohastic), SignalKind::Stohastic);
                                auto dx1_new = 0.25 / inst_freq_first.interpolate(idx_new, SignalKind::Stohastic);

                                auto dx2_new = 0.25 / inst_freq_second.interpolate(idx_new, SignalKind::Stohastic);
                                auto dy2_new = dy1_new;
                                return data.interpolate(idx_new, SignalKind::Stohastic) + dy2_new / dx1_new * dx2_new ;
                            }
                        }

                        return (data.interpolate(idx + 0.25 / inst_freq_second.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic) +
                                data.interpolate(idx - 0.25 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) / 2.0;
                    };
                    for (auto i = 0; i < data.size(); i++) {
                        out[i] = val_expression(i);
                    }
                } else {
                    std::unreachable();
                }
            }
            else if constexpr (filtering_type_k == FilteringType::AverageBased) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = 0;
                        for (double j = -0.5 / (*inst_freq)[i]; j <= 0.5 / (*inst_freq)[i]; j = j + step / (*inst_freq)[i]) {
                            out[i] += data.interpolate(i + j, SignalKind::Stohastic) * step;
                        }
                    }
                } else if (inst_freq_kind == InstFreqKind::Double) {
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = 0;
                        //auto step_left = step*(1./inst_freq[i].left + 1./inst_freq[i].right)/1./inst_freq[i].left;
                        for (double j = -0.5 / (*inst_freq)[i].first; j < 0; j = j + step / (*inst_freq[i]).first) {
                            out[i] += data.interpolate(i + j, SignalKind::Stohastic) * step;
                        }
                        for (double j = 0.0; j < 0.5 / (*inst_freq[i]).second; j = j + step / (*inst_freq[i]).second) {
                            out[i] += data.interpolate(i + j, SignalKind::Stohastic) * step;
                        }
                        //todo left and right steps
                    }
                }
            } else if constexpr (filtering_type_k == FilteringType::Median) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    std::vector<T> buffer;
                    for (int i = 0; i < data.size(); i++) {
                        buffer.clear();
                        for (double j = -0.5 / (*inst_freq)[i]; j <= 0.5 / (*inst_freq)[i]; j = j + step / (*inst_freq)[i]) {
                            buffer.push_back(data.interpolate(i + j, SignalKind::Stohastic));
                        }
                        std::sort(buffer.begin(), buffer.end());
                        if (buffer.size() % 2 == 0) {
                            out[i] = buffer[buffer.size() / 2];
                        } else {
                            out[i] = (buffer[buffer.size() / 2] + buffer[buffer.size() / 2 + 1]) / 2.0;
                        }
                    }
                } else if (inst_freq_kind == InstFreqKind::Double) {
                    std::vector<T> buffer;
                    for (int i = 0; i < data.size(); i++) {
                        buffer.clear();
                        //auto step_left = step*(1./inst_freq[i].left + 1./inst_freq[i].right)/1./inst_freq[i].left;
                        for (double j = -0.5 / (*inst_freq)[i].first; j < 0; j = j + step / (*inst_freq)[i].first) {
                            buffer.push_back(data.interpolate(i + j, SignalKind::Stohastic));
                        }
                        for (double j = 0.0; j < 0.5 / (*inst_freq)[i].second; j = j + step / (*inst_freq)[i].second) {
                            buffer.push_back(data.interpolate(i + j, SignalKind::Stohastic));
                        }
                        std::sort(buffer.begin(), buffer.end());
                        if (buffer.size() % 2 == 0) {
                            out[i] = buffer[buffer.size() / 2];
                        } else {
                            out[i] = (buffer[buffer.size() / 2] + buffer[buffer.size() / 2 + 1]) / 2.0;
                        }
                        //todo left and right steps
                    }
                }
            } else {
                std::unreachable();
            }
        }
    };

    export
    enum class PhaseComputingKind { extremums_based_non_opt, arctg_scaled };

    export
    template<typename U, Filter<U> FilterT, InstFreqComputer<U> InstFreqComputerT,
        PhaseComputer<U> PhaseComputerT, InstFreqComputer<U> InstFreqComputerForModeT, 
        PhaseComputer<U> PhaseComputerForModeT>
    //INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind inst_freq_k>
    struct OptPeriodBasedFilter {
        using AdditionalDataType = SignalPrototype<U>;

        //SimpleVecWrapper<T> mode_base;
        //SimpleVecWrapper<T> inst_freq_buffer2_base;

        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;
        BufferT mode;
        BufferT inst_freq_buffer2;

        double error_threshold = 0.001;
        size_t iter_number = 0;
        size_t good_iter_number = 0;
        size_t true_iter_number = 0;

        std::vector<double> inst_freq_cache;
        double error_old;

        //constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
        constexpr static bool is_filter = true;

        FilterT * filter;
        InstFreqComputerT * inst_freq_computer;
        PhaseComputerT * phase_computer;
        InstFreqComputerForModeT * inst_freq_computer_for_mode;
        PhaseComputerForModeT * phase_computer_for_mode;

        OptPeriodBasedFilter(FilterT & filter,
            InstFreqComputerT & inst_freq_computer,
            PhaseComputerT & phase_computer, 
            InstFreqComputerForModeT & inst_freq_computer_for_mode, 
            PhaseComputerForModeT & phase_computer_for_mode) {
            this->filter = &filter;
            this->inst_freq_computer = &inst_freq_computer;
            this->phase_computer = &phase_computer;
            this->inst_freq_computer_for_mode = &inst_freq_computer_for_mode;
            this->phase_computer_for_mode = &phase_computer_for_mode;
        }

        ~OptPeriodBasedFilter() {
            //delete mode;
            //delete inst_freq_buffer2;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        bool computeIter(const DataType& data, OutType& out, InstFreqType& inst_freq_buffer) {
            using T = typename OutType::SampleType;
            //IC(iter_number);
            iter_number++;
            true_iter_number++;
            filter->compute(data, out, &inst_freq_buffer);

            //after computing filtering
            if (mode.size() != data.size()) {
                static_cast<SimpleVecWrapper<T> *>(mode.base)->vec->clear();
                for (auto i = 0; i < data.size(); i++) {
                    static_cast<SimpleVecWrapper<T> *>(mode.base)->vec->push_back(data[i] - out[i]);
                }
            } else {
                for (auto i = 0; i < data.size(); i++) {
                    mode[i] = data[i] - out[i];
                }
            }
            if (inst_freq_buffer2.size() != inst_freq_buffer.size()) {
                static_cast<SimpleVecWrapper<T> *>(inst_freq_buffer2.base)->vec->clear();
                for (int i = 0; i < inst_freq_buffer.size(); i++) {
                    static_cast<SimpleVecWrapper<T> *>(inst_freq_buffer2.base)->vec->push_back(0.0);
                }
            }

            if (inst_freq_computer_for_mode->is_phase_based()) {
                //todo
                std::unreachable();
            }
            else {
                if constexpr (std::is_same_v<typename InstFreqComputerForModeT::AdditionalDataType, GENERAL::Nil>) {
                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, nullptr);
                } else {
                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, &out);
                }
            }

            double error = UTILITY_MATH::signalsL2NormedDistance
                    <double, InstFreqType, decltype(inst_freq_buffer2)>
                    (inst_freq_buffer, inst_freq_buffer2);

            //double error_old_new = UTILITY_MATH::signalsL2Distance
            //<double, InstFreqType, decltype(inst_freq_cache)>
            //(inst_freq_buffer, inst_freq_cache);
            //IC(error, error_old);
            if (error_old < error) {
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 400) {
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = data[i] - mode[i];
                    }
                    return false;
                }
                if (iter_number - good_iter_number > 4) {
                    for (int i = 0; i < inst_freq_buffer.size(); i++) {
                        inst_freq_buffer[i] = inst_freq_cache[i];
                    }

                    iter_number = good_iter_number;
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = data[i] - mode[i];
                    }
                    //IC(error);
                    auto res = computeIter(data, out, inst_freq_buffer);

                    return res;
                    error = error_old;
                }
            
            } else {
                //IC(error, good_iter_number, iter_number, true_iter_number);
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 4) {
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = data[i] - mode[i];
                    }
                    return false;
                }
                
                good_iter_number = iter_number;
                error_old = error;
            }

            if (error < error_threshold) {
                for (int i = 0; i < data.size(); i++) {
                    out[i] = data[i] - mode[i];
                }
                return false;
            }
            if (iter_number == good_iter_number) {
                for (int i = 0; i < inst_freq_buffer.size(); i++) {
                    inst_freq_cache[i] = inst_freq_buffer[i];
                }
            }

            for (int i = 0; i < inst_freq_buffer.size(); i += std::rand() % (inst_freq_buffer.size() / 5)) {
                //= i + (std::rand() + 1) % (inst_freq_buffer.size()/5)){
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                //double coeff1 = 1.0;
                //double coeff2 = 1.0;
                inst_freq_buffer[i] =
                        (inst_freq_buffer2[i] * coeff2 + inst_freq_buffer[i] * coeff1)
                        / (coeff1 + coeff2);
            }

            for (int i = 0; i < data.size(); i++) {
                out[i] = data[i] - mode[i];
            }

            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        void compute(const DataType& data, OutType& out, InstFreqType * inst_freq_buffer) {
            using T = typename OutType::SampleType;
            //inst_freq_computer->compute(data, inst_freq_buffer, out);
            if (inst_freq_computer->is_phase_based()) {
                //todo
                std::unreachable();
            }
            if constexpr (std::is_same_v<typename InstFreqComputerT::AdditionalDataType, GENERAL::Nil>) {
                inst_freq_computer->compute(data, *inst_freq_buffer, nullptr);
            } else {
                inst_freq_computer->compute(data, *inst_freq_buffer, &out);
            }

            filter->compute(data, out, inst_freq_buffer);

            //after computing filtering
            if (mode.size() != data.size()) {
                static_cast<SimpleVecWrapper<T> *>(mode.base)->vec->clear();
                for (auto i = 0; i < data.size(); i++) {
                    static_cast<SimpleVecWrapper<T> *>(mode.base)->vec->push_back(data[i] - out[i]);
                }
            } else {
                for (auto i = 0; i < data.size(); i++) {
                    mode[i] = data[i] - out[i];
                }
            }
            if (inst_freq_buffer2.size() != inst_freq_buffer->size()) {
                static_cast<SimpleVecWrapper<T> *>(inst_freq_buffer2.base)->vec->clear();
                for (int i = 0; i < inst_freq_buffer->size(); i++) {
                    static_cast<SimpleVecWrapper<T> *>(inst_freq_buffer2.base)->vec->push_back(0.0);
                }
            }


            //inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, out);
            if constexpr (std::is_same_v<typename InstFreqComputerForModeT::AdditionalDataType, GENERAL::Nil>) {
                inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, nullptr);
            } else {
                inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, &out);
            }

            double error = UTILITY_MATH::signalsL2NormedDistance
                    <double, InstFreqType, decltype(inst_freq_buffer2)>
                    (*inst_freq_buffer, inst_freq_buffer2);
            error_old = error;
            //IC(error);

            if (error < error_threshold) {
                for (int i = 0; i < data.size(); i++) {
                    out[i] = data[i] - mode[i];
                }
                return;
            }
            inst_freq_cache.clear();
            for (int i = 0; i < inst_freq_buffer->size(); i++) {
                inst_freq_cache.push_back((*inst_freq_buffer)[i]);
            }
            for (int i = 0; i < inst_freq_buffer->size(); i += std::rand() % (inst_freq_buffer->size() / 5)) {
                //= i + (std::rand() + 1) % (inst_freq_buffer.size()/5)){
                //(*inst_freq_buffer)[i] =
                //        (inst_freq_buffer2[i] + (*inst_freq_buffer)[i] * 10.) / 11.0;
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                //double coeff1 = 1.0;
                //double coeff2 = 1.0;
                (*inst_freq_buffer)[i] =
                        (inst_freq_buffer2[i] * coeff2 + (*inst_freq_buffer)[i] * coeff1)
                        / (coeff1 + coeff2);
            }
            while (computeIter(data, out, *inst_freq_buffer)) {
            }
        }
    };


    export
    template<typename U, Filter<U> FilterT, InstFreqComputer<U> InstFreqComputerT,
        PhaseComputer<U> PhaseComputerT, InstFreqComputer<U> InstFreqComputerForModeT, 
        PhaseComputer<U> PhaseComputerForModeT>
    struct OptPeriodBasedFilterInstFreqDouble {
        using AdditionalDataType = SignalPrototype<U>;


        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;
        BufferT mode;
        GenericSignal<SimpleVecWrapper<std::pair<U, U>>, true> inst_freq_buffer2;

        double error_threshold = 0.1;
        size_t iter_number = 0;
        size_t good_iter_number = 0;
        size_t true_iter_number = 0;

        std::vector<std::pair<double, double>> inst_freq_cache;
        double error_old;

        constexpr static bool is_filter = true;

        FilterT * filter;
        InstFreqComputerT * inst_freq_computer;
        PhaseComputerT * phase_computer;
        InstFreqComputerForModeT * inst_freq_computer_for_mode;
        PhaseComputerForModeT * phase_computer_for_mode;

        OptPeriodBasedFilterInstFreqDouble(FilterT & filter,
            InstFreqComputerT & inst_freq_computer,
            PhaseComputerT & phase_computer, 
            InstFreqComputerForModeT & inst_freq_computer_for_mode, 
            PhaseComputerForModeT & phase_computer_for_mode) {
            this->filter = &filter;
            this->inst_freq_computer = &inst_freq_computer;
            this->phase_computer = &phase_computer;
            this->inst_freq_computer_for_mode = &inst_freq_computer_for_mode;
            this->phase_computer_for_mode = &phase_computer_for_mode;
        }

        ~OptPeriodBasedFilterInstFreqDouble() {
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        bool computeIter(const DataType& data, OutType& out, InstFreqType& inst_freq_buffer) {
            using T = typename OutType::SampleType;
            iter_number++;
            true_iter_number++;
            filter->compute(data, out, &inst_freq_buffer);

            if (mode.size() != data.size()) {
                (mode.base)->vec->clear();
                for (auto i = 0; i < data.size(); i++) {
                    (mode.base)->vec->push_back(data[i] - out[i]);
                }
            } else {
                for (auto i = 0; i < data.size(); i++) {
                    mode[i] = data[i] - out[i];
                }
            }
            if (inst_freq_buffer2.size() != inst_freq_buffer.size()) {
                (inst_freq_buffer2.base)->vec->clear();
                for (int i = 0; i < inst_freq_buffer.size(); i++) {
                    (inst_freq_buffer2.base)->vec->push_back({0.0, 0.0});
                }
            }

            if (inst_freq_computer_for_mode->is_phase_based()) {
                std::unreachable();
            }
            else {
                if constexpr (std::is_same_v<typename InstFreqComputerForModeT::AdditionalDataType, GENERAL::Nil>) {
                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, nullptr);
                } else {
                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, &out);
                }
            }

            double error = UTILITY_MATH::signalsL2NormedDistanceDouble
                    <double, InstFreqType, decltype(inst_freq_buffer2)>
                    (inst_freq_buffer, inst_freq_buffer2);

            if (error_old < error) {
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 400) {
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = data[i] - mode[i];
                    }
                    return false;
                }
                if (iter_number - good_iter_number > 4) {
                    for (int i = 0; i < inst_freq_buffer.size(); i++) {
                        inst_freq_buffer[i] = inst_freq_cache[i];
                    }

                    iter_number = good_iter_number;
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = data[i] - mode[i];
                    }
                    auto res = computeIter(data, out, inst_freq_buffer);

                    return res;
                    error = error_old;
                }
            
            } else {
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 4) {
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = data[i] - mode[i];
                    }
                    return false;
                }
                
                good_iter_number = iter_number;
                error_old = error;
            }

            if (error < error_threshold) {
                for (int i = 0; i < data.size(); i++) {
                    out[i] = data[i] - mode[i];
                }
                return false;
            }

            if (iter_number == good_iter_number) {
                inst_freq_cache.clear();
                for (int i = 0; i < inst_freq_buffer.size(); i++) {
                    inst_freq_cache.push_back(inst_freq_buffer[i]);
                }
            }

            for (int i = 0; i < inst_freq_buffer.size(); i += std::rand() % (inst_freq_buffer.size() / 5)) {
                double coeff = (std::rand() % 500 + 1) / 10.0;
                inst_freq_buffer[i].first =
                        (inst_freq_buffer2[i].first + inst_freq_buffer[i].first * coeff)
                        / (coeff + 1.0);
                inst_freq_buffer[i].second =
                        (inst_freq_buffer2[i].second + inst_freq_buffer[i].second * coeff)
                        / (coeff + 1.0);
            }

            for (int i = 0; i < data.size(); i++) {
                out[i] = data[i] - mode[i];
            }

            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        void compute(const DataType& data, OutType& out, InstFreqType * inst_freq_buffer) {
            using T = typename OutType::SampleType;
            if (inst_freq_computer->is_phase_based()) {
                std::unreachable();
            }
            if constexpr (std::is_same_v<typename InstFreqComputerT::AdditionalDataType, GENERAL::Nil>) {
                inst_freq_computer->compute(data, *inst_freq_buffer, nullptr);
            } else {
                inst_freq_computer->compute(data, *inst_freq_buffer, &out);
            }

            filter->compute(data, out, inst_freq_buffer);

            if (mode.size() != data.size()) {
                (mode.base)->vec->clear();
                for (auto i = 0; i < data.size(); i++) {
                    (mode.base)->vec->push_back(data[i] - out[i]);
                }
            } else {
                for (auto i = 0; i < data.size(); i++) {
                    mode[i] = data[i] - out[i];
                }
            }
            if (inst_freq_buffer2.size() != inst_freq_buffer->size()) {
                (inst_freq_buffer2.base)->vec->clear();
                for (int i = 0; i < inst_freq_buffer->size(); i++) {
                    (inst_freq_buffer2.base)->vec->push_back({0.0, 0.0});
                }
            }

            if constexpr (std::is_same_v<typename InstFreqComputerForModeT::AdditionalDataType, GENERAL::Nil>) {
                inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, nullptr);
            } else {
                inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, &out);
            }

            double error = UTILITY_MATH::signalsL2NormedDistanceDouble
                    <double, InstFreqType, decltype(inst_freq_buffer2)>
                    (*inst_freq_buffer, inst_freq_buffer2);
            error_old = error;
            if (error < error_threshold) {
                for (int i = 0; i < data.size(); i++) {
                    out[i] = data[i] - mode[i];
                }
                return;
            }
            inst_freq_cache.clear();
            for (int i = 0; i < inst_freq_buffer->size(); i++) {
                inst_freq_cache.push_back((*inst_freq_buffer)[i]);
            }
            for (int i = 0; i < inst_freq_buffer->size(); i += std::rand() % (inst_freq_buffer->size() / 5)) {
                (*inst_freq_buffer)[i].first =
                        (inst_freq_buffer2[i].first  + (*inst_freq_buffer)[i].first  * 10.) / 11.0;
                (*inst_freq_buffer)[i].second =
                    (inst_freq_buffer2[i].second  + (*inst_freq_buffer)[i].second  * 10.) / 11.0;
            }
            while (computeIter(data, out, *inst_freq_buffer)) {
            }
        }
    };


    export
    template<typename U, Integrator<U> IntegratorT, Derivator<U> DerivatorT, Filter<U> FilterT, 
        InstFreqComputer<U> InstFreqComputerT, 
            PhaseComputer<U> PhaseComputerT, InstAmplComputer<U> InstAmplComputerT,
                InstFreqComputer<U> InstFreqComputerForModeT, 
                    PhaseComputer<U> PhaseComputerForModeT>
    struct RecursiveFilterInstAmplChanges{
        using AdditionalDataType = SignalPrototype<U>;

        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;

        BufferT prediction_mode;
        BufferT prediction_mode_inst_freq;
        BufferT prediction_mode_ampl;
        BufferT new_result;

        static BufferT & prediction_phase(){
            static BufferT * prediction_mode_phase = new BufferT;
            return *prediction_mode_phase;
        }

        //BufferT computer_buffer;

        std::vector<double> result_cache;
        
        constexpr static bool is_filter = true;

        double error;
        double error_treshhold = 0.1;
        double error_treshhold_muller = 0.1;

        size_t good_iter_number = 0;
        size_t true_iter_number = 0;
        size_t iter_number = 0;

        double error_old = 0.0;


        FilterT * filter;
        InstFreqComputerT * inst_freq_computer;
        InstFreqComputerForModeT * inst_freq_computer_for_mode;
        PhaseComputerT * phase_computer;
        PhaseComputerForModeT * phase_computer_for_mode;
        InstAmplComputerT * inst_ampl_computer;
        INST_AMPL_COMPUTERS::InstAmplNormalizator<double, IntegratorT, DerivatorT, InstAmplComputerT> 
            * inst_ampl_normalizer;

        RecursiveFilterInstAmplChanges(IntegratorT & integrator, DerivatorT & derivator, 
            FilterT & filter, InstFreqComputerT & inst_freq_computer,
                PhaseComputerT & phase_computer, InstAmplComputerT & inst_ampl_computer, 
                    InstFreqComputerForModeT & inst_freq_computer_for_mode,
                        PhaseComputerForModeT & phase_computer_for_mode){
            this->filter = &filter;
            this->inst_freq_computer = &inst_freq_computer;
            this->inst_freq_computer_for_mode = &inst_freq_computer_for_mode;
            this->phase_computer = &phase_computer;
            this->phase_computer_for_mode = &phase_computer_for_mode;
            inst_ampl_normalizer = new INST_AMPL_COMPUTERS::InstAmplNormalizator
                <double, IntegratorT, DerivatorT, InstAmplComputerT> (integrator, derivator, inst_ampl_computer);
        }

        ~RecursiveFilterInstAmplChanges(){
            delete inst_ampl_normalizer;
        }

        template<Signal DataType, Signal OutType, Signal OldResultType>
        bool computeIter(const DataType& data, OutType& result_buffer, OldResultType* computer_buffer){
            iter_number++;
            true_iter_number++;
            using T = typename OutType::SampleType;

            if(result_cache.size()!=data.size()){
                result_cache.clear();
                for(auto i = 0; i < data.size(); i++){
                    result_cache.push_back(0.0);
                }
            }
            
            if(prediction_mode.size() != data.size()){
                prediction_mode.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_inst_freq.size() != data.size()){
                prediction_mode_inst_freq.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_inst_freq.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_ampl.size() != data.size()){
                prediction_mode_ampl.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_ampl.base->vec->push_back(0.0);
                }
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_inst_freq[i] = 0.0;
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_ampl[i] = 0.0;
            }
            for (auto i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = 0.0;
            }
            for(auto i = 0; i < data.size(); i++){
                prediction_mode[i] = data[i] - result_buffer[i];
            }

            if constexpr (InstFreqComputerForModeT::is_phase_based()){
                if (prediction_phase().size()!= data.size){
                    prediction_phase().base->vec->clear();
                    for (auto i = 0; i < data.size(); i++){
                        prediction_phase().base->vec->push_back(0.0);
                    }
                }
                for (auto i = 0; i < data.size(); i++){
                    prediction_phase()[i] = 0.0;
                }
                phase_computer->compute(prediction_mode, prediction_phase(), computer_buffer);
                inst_freq_computer_for_mode->compute(prediction_phase(), prediction_mode_inst_freq, computer_buffer);
            }
            else{
                inst_freq_computer_for_mode->compute(prediction_mode, prediction_mode_inst_freq, computer_buffer);
            }
            
            
            inst_ampl_normalizer->compute(prediction_mode, *computer_buffer, prediction_mode_ampl);

            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = result_buffer[i] + (*computer_buffer)[i];
            }

            filter->compute(*computer_buffer, prediction_mode, &prediction_mode_inst_freq);
            error = UTILITY_MATH::signalsL2NormedDistance
                    <double, decltype(result_buffer), decltype(prediction_mode)>
                    (result_buffer, prediction_mode);

            if (error_old <= error){
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    return false;
                }
                if(iter_number - good_iter_number > 9){
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    iter_number = good_iter_number;
                    return computeIter(data, result_buffer, computer_buffer);
                }
            }
            else if(error < error_old){
                //IC(error_old, error, error_treshhold);
                error_old = error;
                good_iter_number = iter_number;
                
                for (int i = 0; i < data.size(); i++){
                    result_cache[i] = result_buffer[i];
                }
                if(static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    return false;
                }

            }

            
            

            int coeff1 = std::rand() % 30 + 1;
            int coeff2 = std::rand() % 30 + 1;
            for(int i = 0; i < data.size(); i++){
                result_buffer[i] = (result_buffer[i] * coeff1 + prediction_mode[i] * coeff2) / (coeff1 + coeff2);
            }
            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqT>
        void compute(const DataType& data, OutType& result_buffer, InstFreqT* inst_freq){
            iter_number = 0;
            good_iter_number = 0;
            true_iter_number = 0;
            
            //this->inst_freq_computer->compute(data, *inst_freq, &result_buffer);
            this->filter->compute(data, result_buffer, inst_freq);

            double der_avg = 0.0;
            auto sum = 0.0;
            for (int i = 1; i < data.size(); i++){
                der_avg += std::abs((data[i] - result_buffer[i])-(data[i-1] - result_buffer[i-1])) / data.size();
                sum += std::abs(data[i]);
            }
            error_old = sum * data.size();
            der_avg = der_avg / data.size();

            error_treshhold = der_avg * error_treshhold_muller;
            
            computeIter(data, result_buffer, inst_freq);

            while(error > error_treshhold){
                if(!computeIter(data, result_buffer, inst_freq)) break;
            }
        }   
    };



    export
    template<typename U, Filter<U> FilterFirstT, 
        Filter<U> FilterSeconT>
    struct CascadeFilter{
        using AdditionalDataType = SignalPrototype<double>;
        using IdxType = size_t;
        constexpr static bool is_filter = true;

        FilterFirstT * filter_first;
        FilterSeconT * filter_second;

        CascadeFilter(FilterFirstT & filter_first,
            FilterSeconT & filter_second){

            this->filter_first = &filter_first;
            this->filter_second = &filter_second;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        void compute(const DataType& data, OutType& out, 
            const InstFreqType * inst_freq){
            using T = typename OutType::SampleType;
            
            filter_first->compute(data, out, inst_freq);
            auto mode_val = [&](size_t idx){
                return data[idx] - out[idx];
            };
            auto mode_size = [&](){
                return data.size();
            };
            ExpressionWrapper<T, size_t, 
                decltype(mode_val), GENERAL::Nil, decltype(mode_size), false> 
                    mode_expr(mode_val, mode_size);
            GenericSignal<decltype(mode_expr), false> mode(mode_expr);

            GenericSignal<SimpleVecWrapper<T>, true> buffer;

            for(size_t i = 0 ; i < data.size(); i++){
                buffer.base->vec->push_back(0.);
            }

            filter_second->compute(mode, buffer, inst_freq);

            for (int i = 0; i < data.size(); i++){
                out[i] = data[i] - (mode[i] - buffer[i]);
            }
        }
    };

    export
    template<typename U, Integrator<U> IntegratorT, Derivator<U> DerivatorT, Filter<U> FilterT, 
        InstAmplComputer<U> InstAmplComputerT>
    struct RecursiveFilterInstAmplChangesFithConstInstFreq{
        using AdditionalDataType = SignalPrototype<U>;

        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;

        BufferT prediction_mode;
        BufferT prediction_mode_inst_freq;
        BufferT prediction_mode_ampl;
        BufferT new_result;

        static BufferT & prediction_phase(){
            static BufferT * prediction_mode_phase = new BufferT;
            return *prediction_mode_phase;
        }

        //BufferT computer_buffer;

        std::vector<double> result_cache;
        
        constexpr static bool is_filter = true;

        double error;
        double error_treshhold = 0.1;
        double error_treshhold_muller = 0.1;

        size_t good_iter_number = 0;
        size_t true_iter_number = 0;
        size_t iter_number = 0;

        double error_old = 0.0;


        FilterT * filter;
        InstAmplComputerT * inst_ampl_computer;
        INST_AMPL_COMPUTERS::InstAmplNormalizator<double, IntegratorT, DerivatorT, InstAmplComputerT> 
            * inst_ampl_normalizer;

        RecursiveFilterInstAmplChangesFithConstInstFreq(IntegratorT & integrator, DerivatorT & derivator, 
            FilterT & filter, InstAmplComputerT & inst_ampl_computer){
            this->filter = &filter;
            inst_ampl_normalizer = new INST_AMPL_COMPUTERS::InstAmplNormalizator
                <double, IntegratorT, DerivatorT, InstAmplComputerT> (integrator, derivator, inst_ampl_computer);
        }

        ~RecursiveFilterInstAmplChangesFithConstInstFreq(){
            delete inst_ampl_normalizer;
        }

        template<Signal DataType, Signal OutType, Signal OldResultType>
        bool computeIter(const DataType& data, OutType& result_buffer, OldResultType* computer_buffer){
            iter_number++;
            true_iter_number++;
            using T = typename OutType::SampleType;

            if(result_cache.size()!=data.size()){
                result_cache.clear();
                for(auto i = 0; i < data.size(); i++){
                    result_cache.push_back(0.0);
                }
            }
            
            if(prediction_mode.size() != data.size()){
                prediction_mode.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_inst_freq.size() != data.size()){
                prediction_mode_inst_freq.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_inst_freq.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_ampl.size() != data.size()){
                prediction_mode_ampl.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_ampl.base->vec->push_back(0.0);
                }
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_inst_freq[i] = 0.0;
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_ampl[i] = 0.0;
            }
            for(auto i = 0; i < data.size(); i++){
                prediction_mode[i] = data[i] - result_buffer[i];
            }

            for(int i = 0; i < data.size(); i++){
                prediction_mode_inst_freq[i] = (*computer_buffer)[i];
            }
            //auto label = 0;
            //IC(++label);
            inst_ampl_normalizer->compute(prediction_mode, *computer_buffer, prediction_mode_ampl);

            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = result_buffer[i] + (*computer_buffer)[i];
            }
            //IC(++label);
            filter->compute(*computer_buffer, prediction_mode, &prediction_mode_inst_freq);
            //IC(++label);
            error = UTILITY_MATH::signalsL2NormedDistance
                    <double, decltype(result_buffer), decltype(prediction_mode)>
                    (result_buffer, prediction_mode);
            //IC(error);

            if (error_old <= error){
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    for(int i = 0; i < data.size(); i++){
                        (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return false;
                }
                if(iter_number - good_iter_number > 9){
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    iter_number = good_iter_number;
                    for(int i = 0; i < data.size(); i++){
                        (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return computeIter(data, result_buffer, computer_buffer);
                }
            }
            else if(error < error_old){
                //IC(error_old, error, error_treshhold);
                error_old = error;
                good_iter_number = iter_number;
                
                for (int i = 0; i < data.size(); i++){
                    result_cache[i] = result_buffer[i];
                }
                if(static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return false;
                }

            }

            
            

            int coeff1 = std::rand() % 30 + 1;
            int coeff2 = std::rand() % 30 + 1;
            for(int i = 0; i < data.size(); i++){
                result_buffer[i] = (result_buffer[i] * coeff1 + prediction_mode[i] * coeff2) / (coeff1 + coeff2);
            }
            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = prediction_mode_inst_freq[i];
            }
            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqT>
        void compute(const DataType& data, OutType& result_buffer, InstFreqT* inst_freq){
            iter_number = 0;
            good_iter_number = 0;
            true_iter_number = 0;
            
            //this->inst_freq_computer->compute(data, *inst_freq, &result_buffer);
            this->filter->compute(data, result_buffer, inst_freq);

            double der_avg = 0.0;
            auto sum = 0.0;
            for (int i = 1; i < data.size(); i++){
                der_avg += std::abs((data[i] - result_buffer[i])-(data[i-1] - result_buffer[i-1])) / data.size();
                sum += std::abs(data[i]);
            }
            error_old = sum * data.size();
            der_avg = der_avg / data.size();

            error_treshhold = der_avg * error_treshhold_muller;
            
            computeIter(data, result_buffer, inst_freq);

            while(error > error_treshhold){
                if(!computeIter(data, result_buffer, inst_freq)) break;
            }
        }   
    };

    export
    template<typename U, Integrator<U> IntegratorT, Derivator<U> DerivatorT, Filter<U> FilterT, 
        InstAmplComputer<U> InstAmplComputerT>
    struct RecursiveFilterInstAmplChangesFithConstInstFreqDouble{
        //todo
        using AdditionalDataType = SignalPrototype<U>;

        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;

        BufferT prediction_mode;
        BufferT prediction_mode_inst_freq;
        BufferT prediction_mode_ampl;
        BufferT new_result;

        static BufferT & prediction_phase(){
            static BufferT * prediction_mode_phase = new BufferT;
            return *prediction_mode_phase;
        }

        //BufferT computer_buffer;

        std::vector<double> result_cache;
        
        
        constexpr static bool is_filter = true;

        double error;
        double error_treshhold = 0.1;
        double error_treshhold_muller = 0.1;

        size_t good_iter_number = 0;
        size_t true_iter_number = 0;
        size_t iter_number = 0;

        double error_old = 0.0;


        FilterT * filter;
        InstAmplComputerT * inst_ampl_computer;
        INST_AMPL_COMPUTERS::InstAmplNormalizator<double, IntegratorT, DerivatorT, InstAmplComputerT> 
            * inst_ampl_normalizer;

        RecursiveFilterInstAmplChangesFithConstInstFreqDouble(IntegratorT & integrator, DerivatorT & derivator, 
            FilterT & filter, InstAmplComputerT & inst_ampl_computer){
            this->filter = &filter;
            inst_ampl_normalizer = new INST_AMPL_COMPUTERS::InstAmplNormalizator
                <double, IntegratorT, DerivatorT, InstAmplComputerT> (integrator, derivator, inst_ampl_computer);
        }

        ~RecursiveFilterInstAmplChangesFithConstInstFreqDouble(){
            delete inst_ampl_normalizer;
        }

        template<Signal DataType, Signal OutType, Signal OldResultType>
        bool computeIter(const DataType& data, OutType& result_buffer, OldResultType* computer_buffer){
            iter_number++;
            true_iter_number++;
            using T = typename OutType::SampleType;

            if(result_cache.size()!=data.size()){
                result_cache.clear();
                for(auto i = 0; i < data.size(); i++){
                    result_cache.push_back(0.0);
                }
            }
            
            if(prediction_mode.size() != data.size()){
                prediction_mode.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_inst_freq.size() != data.size()){
                prediction_mode_inst_freq.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_inst_freq.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_ampl.size() != data.size()){
                prediction_mode_ampl.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_ampl.base->vec->push_back(0.0);
                }
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_inst_freq[i] = 0.0;
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_ampl[i] = 0.0;
            }
            for(auto i = 0; i < data.size(); i++){
                prediction_mode[i] = data[i] - result_buffer[i];
            }

            for(int i = 0; i < data.size(); i++){
                prediction_mode_inst_freq[i] = (*computer_buffer)[i];
            }
            //auto label = 0;
            //IC(++label);
            inst_ampl_normalizer->compute(prediction_mode, *computer_buffer, prediction_mode_ampl);

            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = result_buffer[i] + (*computer_buffer)[i];
            }
            //IC(++label);
            filter->compute(*computer_buffer, prediction_mode, &prediction_mode_inst_freq);
            //IC(++label);
            error = UTILITY_MATH::signalsL2NormedDistance
                    <double, decltype(result_buffer), decltype(prediction_mode)>
                    (result_buffer, prediction_mode);
            //IC(error);

            if (error_old <= error){
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    for(int i = 0; i < data.size(); i++){
                        (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return false;
                }
                if(iter_number - good_iter_number > 9){
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    iter_number = good_iter_number;
                    for(int i = 0; i < data.size(); i++){
                        (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return computeIter(data, result_buffer, computer_buffer);
                }
            }
            else if(error < error_old){
                //IC(error_old, error, error_treshhold);
                error_old = error;
                good_iter_number = iter_number;
                
                for (int i = 0; i < data.size(); i++){
                    result_cache[i] = result_buffer[i];
                }
                if(static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return false;
                }

            }

            
            

            int coeff1 = std::rand() % 30 + 1;
            int coeff2 = std::rand() % 30 + 1;
            for(int i = 0; i < data.size(); i++){
                result_buffer[i] = (result_buffer[i] * coeff1 + prediction_mode[i] * coeff2) / (coeff1 + coeff2);
            }
            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = prediction_mode_inst_freq[i];
            }
            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqT>
        void compute(const DataType& data, OutType& result_buffer, InstFreqT* inst_freq){
            iter_number = 0;
            good_iter_number = 0;
            true_iter_number = 0;
            
            //this->inst_freq_computer->compute(data, *inst_freq, &result_buffer);
            this->filter->compute(data, result_buffer, inst_freq);

            double der_avg = 0.0;
            auto sum = 0.0;
            for (int i = 1; i < data.size(); i++){
                der_avg += std::abs((data[i] - result_buffer[i])-(data[i-1] - result_buffer[i-1])) / data.size();
                sum += std::abs(data[i]);
            }
            error_old = sum * data.size();
            der_avg = der_avg / data.size();

            error_treshhold = der_avg * error_treshhold_muller;
            
            computeIter(data, result_buffer, inst_freq);

            while(error > error_treshhold){
                if(!computeIter(data, result_buffer, inst_freq)) break;
            }
        }   
    };
}
