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



namespace NP_DSP{
    namespace ONE_D{
        namespace FILTERS{
            export enum class InstFreqKind{Average, Double}; 
            export enum class FilteringType{DerivativeBased, ValueBased, AverageBased, Median};

            export
            template<Signal DataT, Signal OutT, Signal InstFreqT, FilteringType filtering_type_k,
                    Integrator IntegratorT, InstFreqKind inst_freq_k>
            struct NonOptPeriodBasedFilter{
                using DataType = DataT;
                using OutType = OutT;
                using InstFreqType = InstFreqT;
                using AdditionalDataType = InstFreqType;

                //InstFreqT * inst_freq;
                IntegratorT integrator;
                constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
                constexpr static bool is_filter = true;

                double step = 0.1;

                NonOptPeriodBasedFilter(IntegratorT integrator){
                    //this->inst_freq = &inst_freq;
                    this->integrator = integrator;
                }
                NonOptPeriodBasedFilter(){}

                void compute(const DataType & data, OutType & out, const InstFreqType & inst_freq){
                    if constexpr (filtering_type_k == FilteringType::DerivativeBased){
                        if constexpr (inst_freq_kind == InstFreqKind::Average){
                            auto size_expr = [&](){
                                return data.size();
                            };

                            GENERAL::Nil nil;
                            auto val_expression = [&](typename DataType::IdxType idx){
                                return (data.interpolate(idx + 0.5/inst_freq.interpolate(idx)) - data.interpolate(idx - 0.5/inst_freq.interpolate(idx))) *
                                    inst_freq.interpolate(idx);
                            };

                            ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType, 
                                decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                    expr_wrapper (val_expression, size_expr);

                            GenericSignal<decltype(expr_wrapper), false> expr_signal(expr_wrapper);
                            
                            INTEGRATORS::Riman<decltype(expr_signal), OutType,
                            INTEGRATORS::PolygonType::ByPoint> integrator_new;
                            
                            //todo fix its big crucher
                            integrator_new.compute(expr_signal, out, nil);

                            for (int i = 0; i < data.size(); i++) {
                                out[i] += data[0];
                            }
                        }
                        else if constexpr (inst_freq_kind == InstFreqKind::Double){
                            auto inst_freq_first_val_expression = [&](typename DataType::IdxType idx) {
                                return inst_freq[idx].first;
                            };
                            auto inst_freq_second_val_expression = [&](typename DataType::IdxType idx) {
                                return inst_freq[idx].second;
                            };
                            auto size_expr = [&](){
                                return data.size();
                            };
                            ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_first_val_expression), GENERAL::Nil, decltype(size_expr), false>
                                    inst_freq_first_expr_wrapper (inst_freq_first_val_expression, size_expr);
                            GenericSignal<decltype(inst_freq_first_expr_wrapper), false> inst_freq_first(inst_freq_first_expr_wrapper);

                            ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_second_val_expression), GENERAL::Nil, decltype(size_expr), false>
                                    inst_freq_second_expr_wrapper (inst_freq_second_val_expression, size_expr);
                            GenericSignal<decltype(inst_freq_second_expr_wrapper), false> inst_freq_second(inst_freq_second_expr_wrapper);

                            auto val_expression = [&](typename DataType::IdxType idx) {
                                return (data.interpolate(idx + 0.5/inst_freq_second.interpolate(idx)) -
                                    data.interpolate(idx - 0.5/inst_freq_first.interpolate(idx))) /
                                        (0.5 / inst_freq_second.interpolate(idx) + 0.5 / inst_freq_first.interpolate(idx));
                            };
                            GENERAL::Nil nil;

                            ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                    expr_wrapper (val_expression, size_expr);

                            GenericSignal<decltype(expr_wrapper), false> expr_signal(expr_wrapper);

                            INTEGRATORS::Riman<decltype(expr_signal), OutType,
                            INTEGRATORS::PolygonType::ByPoint> integrator_new;

                            //todo fix its big crucher
                            integrator_new.compute(expr_signal, out, nil);

                            for (int i = 0; i < data.size(); i++) {
                                out[i] += data[0];
                            }
                        }
                        else{
                            std::unreachable();
                        }
                    }
                    else if constexpr (filtering_type_k == FilteringType::ValueBased){
                        if constexpr (inst_freq_kind == InstFreqKind::Average){
                            auto val_expression = [&](typename DataType::IdxType idx){
                                return (data.interpolate(idx + 0.25/inst_freq.interpolate(idx)) 
                                    + data.interpolate(idx - 0.25/inst_freq.interpolate(idx))) / 2;
                            };
                            for (auto i = 0; i < data.size(); i++){
                                out[i] = val_expression(i);
                            }
                        }
                        else if constexpr (inst_freq_kind == InstFreqKind::Double){
                            auto inst_freq_first_val_expression = [&](typename DataType::IdxType idx) {
                                return inst_freq[idx].first;
                            };
                            auto inst_freq_second_val_expression = [&](typename DataType::IdxType idx) {
                                return inst_freq[idx].second;
                            };
                            auto size_expr = [&]() {
                                return data.size();
                            };
                            ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_first_val_expression), GENERAL::Nil, decltype(size_expr), false>
                                    inst_freq_first_expr_wrapper (inst_freq_first_val_expression, size_expr);
                            GenericSignal<decltype(inst_freq_first_expr_wrapper), false> inst_freq_first(inst_freq_first_expr_wrapper);

                            ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_second_val_expression), GENERAL::Nil, decltype(size_expr), false>
                                    inst_freq_second_expr_wrapper (inst_freq_second_val_expression, size_expr);
                            GenericSignal<decltype(inst_freq_second_expr_wrapper), false> inst_freq_second(inst_freq_second_expr_wrapper);

                            auto val_expression = [&](typename DataType::IdxType idx) {
                                return (data.interpolate(idx + 0.25/inst_freq_second.interpolate(idx)) +
                                    data.interpolate(idx - 0.25/inst_freq_first.interpolate(idx))) / 2.0;
                            };
                            for (auto i = 0; i < data.size(); i++){
                                out[i] = val_expression(i);
                            }
                        }
                        else{
                            std::unreachable();
                        }
                    }
                    else if constexpr (filtering_type_k == FilteringType::AverageBased) {
                        if constexpr (inst_freq_kind == InstFreqKind::Average) {
                            for (int i = 0; i < data.size(); i++) {
                                out[i] = 0;
                                for (double j = -0.5/inst_freq[i]; j <= 0.5/inst_freq[i]; j = j + step/inst_freq[i]) {
                                    out[i] += data.interpolate(i + j) * step ;
                                }
                            }
                        }
                        else if (inst_freq_kind == InstFreqKind::Double) {
                            for (int i = 0; i < data.size(); i++) {
                                out[i] = 0;
                                //auto step_left = step*(1./inst_freq[i].left + 1./inst_freq[i].right)/1./inst_freq[i].left;
                                for (double j = -0.5/inst_freq[i].first; j < 0; j = j + step/inst_freq[i].first) {
                                    out[i] += data.interpolate(i + j) * step ;
                                }
                                for (double j = 0.0; j < 0.5/inst_freq[i].second; j = j + step/inst_freq[i].second) {
                                    out[i] += data.interpolate(i + j) * step ;
                                }
                                //todo left and right steps
                            }
                        }
                    }
                    else if constexpr (filtering_type_k == FilteringType::Median) {
                        if constexpr (inst_freq_kind == InstFreqKind::Average) {
                            std::vector<typename DataT::SampleType> buffer;
                            for (int i = 0; i < data.size(); i++) {
                                buffer.clear();
                                for (double j = -0.5/inst_freq[i]; j <= 0.5/inst_freq[i]; j = j + step/inst_freq[i]) {
                                    buffer.push_back(data.interpolate(i + j));
                                }
                                std::sort(buffer.begin(), buffer.end());
                                if (buffer.size()%2 == 0) {
                                    out[i] = buffer[buffer.size()/2];
                                }
                                else {
                                    out[i] = (buffer[buffer.size()/2] + buffer[buffer.size()/2 + 1]) / 2.0;
                                }
                            }

                        }
                        else if (inst_freq_kind == InstFreqKind::Double) {
                            std::vector<typename DataT::SampleType> buffer;
                            for (int i = 0; i < data.size(); i++) {
                                buffer.clear();
                                //auto step_left = step*(1./inst_freq[i].left + 1./inst_freq[i].right)/1./inst_freq[i].left;
                                for (double j = -0.5/inst_freq[i].first; j < 0; j = j + step/inst_freq[i].first) {
                                    buffer.push_back(data.interpolate(i + j));
                                }
                                for (double j = 0.0; j < 0.5/inst_freq[i].second; j = j + step/inst_freq[i].second) {
                                    buffer.push_back(data.interpolate(i + j));
                                }
                                std::sort(buffer.begin(), buffer.end());
                                if (buffer.size()%2 == 0) {
                                    out[i] = buffer[buffer.size()/2];
                                }
                                else {
                                    out[i] = (buffer[buffer.size()/2] + buffer[buffer.size()/2 + 1]) / 2.0;
                                }
                                //todo left and right steps
                            }
                        }
                    }
                    else{
                        std::unreachable();
                    }
                }
            };

            export
            enum class InstFreqComputerKind {extremums_based,
                phase_based_momental, phase_based_time_average,
                phase_based_derive_average};

            export 
            enum class PhaseComputingKind {extremums_based_non_opt, arctg_scaled};

            export
            template<Signal DataT, Signal OutT, Signal InstFreqT, PHASE_COMPUTERS::ExtremumsKind kind_e,
                    FilteringType filtering_type_k,
                    Integrator IntegratorT, Derivator DerivatorT, InstFreqComputerKind inst_freq_computer_k,
                    PhaseComputingKind phase_computer_k>
                    //INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind inst_freq_k>
            struct OptPeriodBasedFilter{
                using DataType = DataT;
                using OutType = OutT;
                using InstFreqType = InstFreqT;
                using AdditionalDataType = InstFreqType;
                IntegratorT integrator;
                DerivatorT derivator;
                GenericSignal<SimpleVecWrapper<double>, true> mode;
                GenericSignal<SimpleVecWrapper<double>, true> inst_freq_buffer2;

                double error_threshold = 0.1;
                using BuffT = decltype(mode);
                size_t iter_number = 0;
                size_t good_iter_number = 0;
                size_t true_iter_number = 0;

                std::vector<double> inst_freq_cache;
                double error_old;

                //constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
                constexpr static bool is_filter = true;
                NonOptPeriodBasedFilter<DataType, OutType,
                        AdditionalDataType, filtering_type_k, 
                            IntegratorT, InstFreqKind::Average> * filter;

                static constexpr bool isSameAddNil(){
                    using PhaseComputerT = PHASE_COMPUTERS::ExtremumsBasedNonOpt
                       <DataT, InstFreqT, kind_e, DerivatorT>;
                    if constexpr (inst_freq_computer_k == InstFreqComputerKind::extremums_based){
                        return std::convertible_to<typename
                            INST_FREQ_COMPUTERS::ExtremumsBased <DataT, InstFreqT,
                                    INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear>::AdditionalDataType
                                    , GENERAL::Nil>;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::phase_based_momental){
                        return std::convertible_to<typename
                        INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT, DerivatorT,
                            INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental, PhaseComputerT>::AdditionalDataType
                                , GENERAL::Nil>;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::phase_based_time_average){
                        return std::convertible_to<typename
                        INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                                PhaseComputerT>::AdditionalDataType
                                    , GENERAL::Nil>;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::phase_based_derive_average){
                        return std::convertible_to<typename
                        INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage,
                                PhaseComputerT>::AdditionalDataType
                                    , GENERAL::Nil>;
                    }
                }

                auto & inst_freq_computer(){
                    if constexpr (inst_freq_computer_k == InstFreqComputerKind::extremums_based){
                        static INST_FREQ_COMPUTERS::ExtremumsBased
                            <DataT, InstFreqT, 
                                INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear>
                                    inst_freq_computer;
                        return inst_freq_computer;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::phase_based_momental){
                        if constexpr (phase_computer_k == PhaseComputingKind::extremums_based_non_opt) {
                            static PHASE_COMPUTERS::ExtremumsBasedNonOpt
                                <DataT, InstFreqT, kind_e, DerivatorT> phase_computer;

                            static INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);

                            return inst_freq_computer;
                        }
                        else if constexpr (phase_computer_k == PhaseComputingKind::arctg_scaled) {
                            static PHASE_COMPUTERS::ArctgScaledToExtremums
                                <DataT, InstFreqT, OutType, kind_e, decltype(integrator),
                                    decltype(derivator)> phase_computer(integrator, derivator);

                            static INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);
                            return inst_freq_computer;
                        }
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::phase_based_time_average){
                        if constexpr (phase_computer_k == PhaseComputingKind::extremums_based_non_opt) {
                            static PHASE_COMPUTERS::ExtremumsBasedNonOpt
                                <DataT, InstFreqT, kind_e, DerivatorT> phase_computer;

                            static INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);

                            return inst_freq_computer;
                        }
                        else if constexpr (phase_computer_k == PhaseComputingKind::arctg_scaled) {
                            static PHASE_COMPUTERS::ArctgScaledToExtremums
                                <DataT, InstFreqT, OutType, kind_e, decltype(integrator),
                                    decltype(derivator)> phase_computer(integrator, derivator);

                            static INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);
                            return inst_freq_computer;
                        }
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::phase_based_derive_average){
                        if constexpr (phase_computer_k == PhaseComputingKind::extremums_based_non_opt) {
                            static PHASE_COMPUTERS::ExtremumsBasedNonOpt
                                <DataT, InstFreqT, kind_e, DerivatorT> phase_computer;

                            static INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);

                            return inst_freq_computer;
                        }
                        else if constexpr (phase_computer_k == PhaseComputingKind::arctg_scaled) {
                            static PHASE_COMPUTERS::ArctgScaledToExtremums
                                <DataT, InstFreqT, OutType, kind_e, decltype(integrator),
                                    decltype(derivator)> phase_computer(integrator, derivator);

                            static INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);
                            return inst_freq_computer;
                        }
                    }
                }

                auto & inst_freq_computer_for_mode(){
                    if constexpr (inst_freq_computer_k == InstFreqComputerKind::extremums_based){
                        static INST_FREQ_COMPUTERS::ExtremumsBased
                            <DataT, InstFreqT,
                                INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear>
                                    inst_freq_computer;
                        return inst_freq_computer;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::phase_based_momental){
                        if constexpr (phase_computer_k == PhaseComputingKind::extremums_based_non_opt) {
                            static PHASE_COMPUTERS::ExtremumsBasedNonOpt
                                <DataT, InstFreqT, kind_e, DerivatorT> phase_computer;

                            static INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);

                            return inst_freq_computer;
                        }
                        else if constexpr (phase_computer_k == PhaseComputingKind::arctg_scaled) {
                            static PHASE_COMPUTERS::ArctgScaledToExtremums
                                <DataT, InstFreqT, OutType, kind_e, decltype(integrator),
                                    decltype(derivator)> phase_computer(integrator, derivator);

                            static INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);
                            return inst_freq_computer;
                        }
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::phase_based_time_average){
                        if constexpr (phase_computer_k == PhaseComputingKind::extremums_based_non_opt) {
                            static PHASE_COMPUTERS::ExtremumsBasedNonOpt
                                <BuffT, BuffT, kind_e, DerivatorT> phase_computer;

                            static INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);

                            return inst_freq_computer;
                        }
                        else if constexpr (phase_computer_k == PhaseComputingKind::arctg_scaled) {
                            static PHASE_COMPUTERS::ArctgScaledToExtremums
                                <BuffT, BuffT, OutType, kind_e, decltype(integrator),
                                    decltype(derivator)> phase_computer(integrator, derivator);

                            static INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);
                            return inst_freq_computer;
                        }
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::phase_based_derive_average){
                        if constexpr (phase_computer_k == PhaseComputingKind::extremums_based_non_opt) {
                            static PHASE_COMPUTERS::ExtremumsBasedNonOpt
                                <BuffT, BuffT, kind_e, DerivatorT> phase_computer;

                            static INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);

                            return inst_freq_computer;
                        }
                        else if constexpr (phase_computer_k == PhaseComputingKind::arctg_scaled) {
                            static PHASE_COMPUTERS::ArctgScaledToExtremums
                                <BuffT, BuffT, OutType, kind_e, decltype(integrator),
                                    decltype(derivator)> phase_computer(integrator, derivator);

                            static INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT,
                                DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage,
                                    decltype(phase_computer)>
                                        inst_freq_computer(integrator, derivator, phase_computer);
                            return inst_freq_computer;
                        }
                    }
                }

                GENERAL::Nil nil;
                
                OptPeriodBasedFilter(IntegratorT integrator, DerivatorT derivator){
                    //this->inst_freq = &inst_freq;
                    this->integrator = integrator;
                    this->derivator = derivator;
                    filter = new NonOptPeriodBasedFilter
                        <DataType, OutType, AdditionalDataType, filtering_type_k, 
                            IntegratorT, InstFreqKind::Average> (integrator);
                }

                ~OptPeriodBasedFilter(){
                    delete filter;
                }

                bool computeIter(const DataType & data, OutType & out, InstFreqType & inst_freq_buffer){
                    iter_number++;
                    true_iter_number++;
                    filter->compute(data, out, inst_freq_buffer);

                    //after computing filtering
                    if (mode.size() != data.size()){
                        mode.base->vec->clear();
                        for (auto i = 0; i < data.size(); i++){
                            mode.base->vec->push_back(data[i] - out[i]);
                        }
                    }
                    else{
                        for (auto i = 0; i < data.size(); i++){
                            mode[i] = data[i] - out[i];
                        }
                    }
                    if(inst_freq_buffer2.size() != inst_freq_buffer.size()){
                        inst_freq_buffer2.base->vec->clear();
                        for (int i = 0; i < inst_freq_buffer.size(); i++){
                            inst_freq_buffer2.base->vec->push_back(0);
                        }
                    }
                    //out.show(PlottingKind::Simple);
                    //mode.show(PlottingKind::Simple);

                    if constexpr (isSameAddNil()){
                        inst_freq_computer_for_mode().compute(mode, inst_freq_buffer2, nil);
                    }
                    else{
                        inst_freq_computer_for_mode().compute(mode, inst_freq_buffer2, out);
                    }

                    double error = UTILITY_MATH::signalsL2Distance
                        <double, InstFreqType, decltype(inst_freq_buffer2)>
                            (inst_freq_buffer, inst_freq_buffer2);

                    //double error_old_new = UTILITY_MATH::signalsL2Distance
                        //<double, InstFreqType, decltype(inst_freq_cache)>
                            //(inst_freq_buffer, inst_freq_cache);
                    //IC(error, error_old);
                    if (error_old < error) {
                        if (iter_number - good_iter_number > 4) {
                            for (int i = 0; i < inst_freq_buffer.size(); i++) {
                                inst_freq_buffer[i] = inst_freq_cache[i];
                            }

                            iter_number = good_iter_number;
                            for (int i = 0; i < data.size(); i++){
                                out[i] = data[i] - mode[i];
                            }
                            auto res = computeIter(data, out, inst_freq_buffer);

                            return res;
                            error = error_old;
                        }
                    }
                    else {
                        if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                            > 4.5) {
                            for (int i = 0; i < data.size(); i++){
                                out[i] = data[i] - mode[i];
                            }
                            return false;
                        }
                        IC(error, good_iter_number, iter_number, true_iter_number);
                        good_iter_number = iter_number;
                        error_old = error;
                    }

                    if (error < error_threshold){
                        for (int i = 0; i < data.size(); i++){
                            out[i] = data[i] - mode[i];
                        }
                        return false;
                    }

                    //inst_freq_buffer.show(PlottingKind::Simple);
                    //inst_freq_buffer2.show(PlottingKind::Simple);

                    if (iter_number == good_iter_number) {
                        inst_freq_cache.clear();
                        for(int i = 0; i < inst_freq_buffer.size(); i++) {
                            inst_freq_cache.push_back(inst_freq_buffer[i]);
                        }
                    }

                    for (int i = 0; i < inst_freq_buffer.size(); i += std::rand() % (inst_freq_buffer.size()/5)){ //= i + std::rand() % (inst_freq_buffer.size()/5)){
                        double coeff = (std::rand() % 500) / 10.0;
                        inst_freq_buffer[i] =
                            (inst_freq_buffer2[i] + inst_freq_buffer[i] * coeff)
                                / (coeff + 1.0);
                        //inst_freq_buffer.show();
                    }

                    for (int i = 0; i < data.size(); i++){
                        out[i] = data[i] - mode[i];
                    }

                    return true;
                }

                void compute(const DataType & data, OutType & out, InstFreqType & inst_freq_buffer) {
                    //inst_freq_computer->compute(data, inst_freq_buffer, out);
                    if constexpr (isSameAddNil()){
                        inst_freq_computer().compute(data, inst_freq_buffer, nil);
                    }
                    else{
                        inst_freq_computer().compute(data, inst_freq_buffer, out);
                    }
                    
                    filter->compute(data, out, inst_freq_buffer);

                    //after computing filtering
                    if (mode.size() != data.size()){
                        mode.base->vec->clear();
                        for (auto i = 0; i < data.size(); i++){
                            mode.base->vec->push_back(data[i] - out[i]);
                        }
                    }
                    else{
                        for (auto i = 0; i < data.size(); i++){
                            mode[i] = data[i] - out[i];
                        }
                    }
                    if(inst_freq_buffer2.size() != inst_freq_buffer.size()){
                        inst_freq_buffer2.base->vec->clear();
                        for (int i = 0; i < inst_freq_buffer.size(); i++){
                            inst_freq_buffer2.base->vec->push_back(0);
                        }
                    }
                    

                    //inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, out);
                    if constexpr (isSameAddNil()){
                        inst_freq_computer_for_mode().compute(mode, inst_freq_buffer2, nil);
                    }
                    else{
                        inst_freq_computer_for_mode().compute(mode, inst_freq_buffer2, out);
                    }

                    double error = UTILITY_MATH::signalsL2Distance
                        <double, InstFreqType, decltype(inst_freq_buffer2)>
                            (inst_freq_buffer, inst_freq_buffer2);
                    error_old = error;
                    //IC(error);
                    
                    //mode.show(PlottingKind::Simple);
                    //inst_freq_buffer.show(PlottingKind::Simple);
                    if (error < error_threshold){
                        for (int i = 0; i < data.size(); i++){
                            out[i] = data[i] - mode[i];
                        }
                        return;
                    }
                    inst_freq_cache.clear();
                    for(int i = 0; i < inst_freq_buffer.size(); i++) {
                        inst_freq_cache.push_back(inst_freq_buffer[i]);
                    }
                    for (int i = 0; i < inst_freq_buffer.size(); i += std::rand() % (inst_freq_buffer.size()/5)){ //= i + std::rand() % (inst_freq_buffer.size()/5)){
                        inst_freq_buffer[i] =
                            (inst_freq_buffer2[i] + inst_freq_buffer[i] * 10.) / 11.0;
                        //inst_freq_buffer.show();
                    }
                    while(computeIter(data, out, inst_freq_buffer)){}
                }
            };
        }
    }
}