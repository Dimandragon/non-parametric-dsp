module;

#include "icecream.hpp"

export module inst_freq_computers;

import <matplot/matplot.h>;
import npdsp_concepts;
import signals;
import derivators;
import integrators;
import <cmath>;
import <numbers>;
import <vector>;
import <utility>;
import utility_math;
import approximators;
import <optional>;
import npdsp_config;
import <complex>;

namespace NP_DSP{
    namespace ONE_D{
        namespace InstFreqComputers{
            export
            enum class InstFreqDerivativeBasedKind {Momental, TimeAverage, DeriveAverage, DeriveDouble};

            export
            template<Signal DataT, Signal OutT, /*SignalWrapper OptFn, OptFn opt_fn,*/
                    Integrator IntegratorT, Derivator DerivatorT, InstFreqDerivativeBasedKind kind>
            struct DerivativeBased{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = OutT;

                constexpr static bool is_inst_freq_computer = true;
                constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

                //using OptFunction = OptFn;
                //OptFunction opt_function = opt_fn;

                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;

                //static_assert(OutType::is_writable == true);

                DerivativeBased(IntegratorT integrator_o,
                                DerivatorT derivator_o)
                {
                    integrator = integrator_o;
                    derivator = derivator_o;
                }

                void compute(const DataType & data, OutType & out, AdditionalDataType & computer_buffer)
                {
                    GENERAL::Nil nil;
                    if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental){
                        derivator.compute(data, computer_buffer, nil);
                        for (int i = 0; i < data.size(); i++){
                            computer_buffer[i] = std::atan(computer_buffer[i]);
                        }
                        derivator.compute(computer_buffer, out, nil);
                        for (int i = 0; i < data.size(); i++){
                            out[i] = std::abs(out[i]) / (std::numbers::pi * 2.0);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage){

                        if constexpr (NP_DSP::CONFIG::debug){
                            //data.show(PlottingKind::Interpolate);
                        }
                        derivator.compute(data, computer_buffer, nil);
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot1";
                            //IC(mark);
                            //computer_buffer.show(PlottingKind::Interpolate);
                        }
                        for (int i = 0; i < data.size(); i++){
                            computer_buffer[i] = std::atan(computer_buffer[i]);
                        }
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot2";
                            //IC(mark);
                            //computer_buffer.show(PlottingKind::Interpolate);
                        }
                        derivator.compute(computer_buffer, out, nil);
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot3";
                            //IC(mark);
                            //out.show(PlottingKind::Interpolate);
                        }
                        for (int i = 0; i < data.size(); i++){
                            out[i] = std::abs(out[i]);
                        }
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot4";
                            //IC(mark);
                            //out.show(PlottingKind::Interpolate);
                        }
                        integrator.compute(out, computer_buffer, nil);
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot5";
                            //IC(mark);
                            //computer_buffer.show(PlottingKind::Interpolate);
                        }
                        
                        for (auto i = 0; i < data.size(); i++){
                            auto approx_answer = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer = approx_answer;
                            auto counter = 0;
                            while (approx_answer < 2.0 * std::numbers::pi){
                                counter++;
                                old_approx_answer = approx_answer;
                                approx_answer = computer_buffer.interpolate(i+counter) - computer_buffer.interpolate(i - counter);
                            }
                            auto left_loss = std::numbers::pi * 2.0 - old_approx_answer;
                            auto right_loss = approx_answer - std::numbers::pi * 2.0;
                            auto sum_loss = left_loss + right_loss;
                            auto period = (static_cast<OutType::SampleType>(counter) - right_loss/sum_loss) * 2;
                            out[i] = static_cast<OutType::SampleType>(1.0/period);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveAverage){
                        derivator.compute(data, computer_buffer, nil);
                        for (int i = 0; i < data.size(); i++){
                            computer_buffer[i] = std::atan(computer_buffer[i]);
                        }
                        derivator.compute(computer_buffer, out, nil);
                        for (int i = 0; i < data.size(); i++){
                            out[i] = std::abs(out[i]);
                        }
                        integrator.compute(out, computer_buffer, nil);
                        if constexpr (NP_DSP::CONFIG::debug){
                            computer_buffer.show(PlottingKind::Interpolate);
                        }
                        //SampleType right_edge;// = computer_buffer.template findMonotone<double>(computer_buffer[0] + std::numbers::pi, {0}, {});
                        //SampleType left_edge;// = computer_buffer.template findMonotone<double>(computer_buffer[0] - std::numbers::pi, {}, {0});
                        //SampleType period; //= right_edge + left_edge;
                        //out[0] = static_cast<OutType::SampleType>(1.0/period);
                        for (auto i = 0; i < data.size(); i++){
                            //right_edge = computer_buffer.template findMonotone<double>(computer_buffer[i] + std::numbers::pi, {right_edge}, {});
                            //left_edge = computer_buffer.template findMonotone<double>(computer_buffer[i] - std::numbers::pi, {left_edge}, {});
                            auto approx_answer_right = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer_right = approx_answer_right;
                            auto counter = 0;
                            while (approx_answer_right < std::numbers::pi){
                                counter++;
                                old_approx_answer_right = approx_answer_right;
                                //approx_answer_right = compute_buffer.interpolate(i+counter) - compute_buffer.interpolate(i);
                                approx_answer_right = computer_buffer.interpolate(i + counter) - computer_buffer.interpolate(i);
                                //std::abs(computer_buffer.interpolate(i+counter)-computer_buffer.interpolate(i+counter-1));
                                //        std::abs(computer_buffer.interpolate(i-counter)-computer_buffer.interpolate(i-counter+1));
                            }
                            auto left_loss_right_edge = std::numbers::pi - old_approx_answer_right;
                            auto right_loss_right_edge = approx_answer_right - std::numbers::pi;
                            auto sum_loss_right_edge = left_loss_right_edge + right_loss_right_edge;
                            auto right_edge = static_cast<OutType::SampleType>(counter) - right_loss_right_edge/sum_loss_right_edge;
                            
                            //auto right_edge = computer_buffer.template findMonotone<double>(computer_buffer[i] + std::numbers::pi, {i + counter-1}, {i + counter});
                            
                            auto approx_answer_left = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer_left = approx_answer_left;
                            counter = 0;
                            while (approx_answer_left < std::numbers::pi){
                                counter++;
                                old_approx_answer_left = approx_answer_left;
                                //approx_answer_right = compute_buffer.interpolate(i+counter) - compute_buffer.interpolate(i);
                                approx_answer_left = computer_buffer.interpolate(i)-computer_buffer.interpolate(i-counter);
                            }
                            //auto left_edge = computer_buffer.template findMonotone<double>(computer_buffer[i] - std::numbers::pi, {i - counter}, {i - counter + 1});
                            auto left_loss_left_edge = std::numbers::pi - old_approx_answer_left;
                            auto right_loss_left_edge = approx_answer_left - std::numbers::pi;
                            auto sum_loss_left_edge = left_loss_left_edge + right_loss_left_edge;
                            auto left_edge = static_cast<OutType::SampleType>(counter) - right_loss_left_edge/sum_loss_left_edge; 

                            auto period = right_edge + left_edge;
                            out[i] = static_cast<OutType::SampleType>(1.0/period);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveDouble){
                        derivator.compute(data, computer_buffer, nil);
                        for (int i = 0; i < data.size(); i++){
                            computer_buffer[i] = std::atan(computer_buffer[i]);
                        }
                        auto out_ref_forward_getter_lamda = [&out](OutType::IdxType idx){
                            auto & ref = out[idx];
                            auto & out_ref = ref.forward;
                            return out_ref;
                        };
                        auto out_ref_backward_getter_lamda = [&out](OutType::IdxType idx){
                            auto & ref = out[idx];
                            auto & out_ref = ref.backward;
                            return out_ref;
                        };
                        auto out_val_getter_lamda = [&out](OutType::IdxType idx){
                            return out[idx].forward;
                        };
                        auto get_size_lamda = [out](){
                            return out.size();
                        };
                        auto signal_base = ExpressionWrapper<typename OutType::SampleType,
                            typename OutType::IdxType, decltype(out_val_getter_lamda), decltype(out_ref_forward_getter_lamda),
                            decltype(get_size_lamda), true>(out_val_getter_lamda, out_ref_forward_getter_lamda, get_size_lamda);
                        auto signal_for_compute = GenericSignal(signal_base);

                        derivator.compute(computer_buffer, signal_for_compute, nil);
                        for (int i = 0; i < data.size(); i++){
                            out_ref_forward_getter_lamda(i) = std::abs(out_val_getter_lamda(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(signal_for_compute, computer_buffer);
                        //if constexpr (NP_DSP::CONFIG::debug){
                        //    computer_buffer.show(PlottingKind::Interpolate);
                        //}
                        //todo

                        auto right_edge = computer_buffer.findIncr(computer_buffer[0] + std::numbers::pi, {0}, nil);
                        auto left_edge = computer_buffer.findIncr(computer_buffer[0] - std::numbers::pi, nil, {0});
                        out_ref_forward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/left_edge);
                        out_ref_backward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/right_edge);
                        for (auto i = 1; i < data.size(); i++){
                            right_edge = computer_buffer.findIncr(computer_buffer[i] + std::numbers::pi, {right_edge}, nil);
                            left_edge = computer_buffer.findIncr(computer_buffer[i] - std::numbers::pi, {left_edge}, nil);
                            out_ref_forward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/left_edge);
                            out_ref_backward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/right_edge);
                        }
                    }
                }
            };

            export
            template<Signal DataT, Signal OutT, /*SignalWrapper OptFn, OptFn opt_fn,*/
                    Integrator IntegratorT, Derivator DerivatorT, InstFreqDerivativeBasedKind kind>
            struct DerivativeBasedWithExternalOptParametr{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = OutT;

                constexpr static bool is_inst_freq_computer = true;
                constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

                //using OptFunction = OptFn;
                //OptFunction opt_function = opt_fn;

                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;

                //static_assert(OutType::is_writable == true);

                DerivativeBasedWithExternalOptParametr(IntegratorT integrator_o,
                                DerivatorT derivator_o)
                {
                    integrator = integrator_o;
                    derivator = derivator_o;
                }

                void compute(const DataType & data, OutType & out, AdditionalDataType & computer_buffer)
                {
                    GENERAL::Nil nil;
                    if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental){
                        derivator.compute(data, out, nil);
                        for (int i = 0; i < data.size(); i++){
                            out[i] = std::atan(out[i]) * computer_buffer[i];
                        }
                        derivator.compute(out, computer_buffer, nil);
                        for (int i = 0; i < data.size(); i++){
                            out[i] = std::abs(computer_buffer[i]) / (std::numbers::pi * 2.0);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage){

                        if constexpr (NP_DSP::CONFIG::debug){
                            //data.show(PlottingKind::Interpolate);
                        }
                        derivator.compute(data, out, nil);
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot1";
                            //IC(mark);
                            //out.show(PlottingKind::Interpolate);
                        }
                        for (int i = 0; i < data.size(); i++){
                            computer_buffer[i] = std::atan(out[i]) * computer_buffer[i];
                        }
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot2";
                            //IC(mark);
                            //computer_buffer.show(PlottingKind::Interpolate);
                        }
                        derivator.compute(computer_buffer, out, nil);
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot3";
                            //IC(mark);
                            //out.show(PlottingKind::Interpolate);
                        }
                        for (int i = 0; i < data.size(); i++){
                            out[i] = std::abs(out[i]);
                        }
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot4";
                            //IC(mark);
                            //out.show(PlottingKind::Interpolate);
                        }
                        integrator.compute(out, computer_buffer, nil);
                        if constexpr (NP_DSP::CONFIG::debug){
                            //std::string mark = "plot5";
                            //IC(mark);
                            //computer_buffer.show(PlottingKind::Interpolate);
                        }
                        
                        for (auto i = 0; i < data.size(); i++){
                            auto approx_answer = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer = approx_answer;
                            auto counter = 0;
                            while (approx_answer < 2.0 * std::numbers::pi){
                                counter++;
                                old_approx_answer = approx_answer;
                                approx_answer = computer_buffer.interpolate(i+counter) - computer_buffer.interpolate(i - counter);
                            }
                            auto left_loss = std::numbers::pi * 2.0 - old_approx_answer;
                            auto right_loss = approx_answer - std::numbers::pi * 2.0;
                            auto sum_loss = left_loss + right_loss;
                            auto period = (static_cast<OutType::SampleType>(counter) - right_loss/sum_loss) * 2;
                            out[i] = static_cast<OutType::SampleType>(1.0/period);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveAverage){
                        derivator.compute(data, out, nil);
                        for (int i = 0; i < data.size(); i++){
                            computer_buffer[i] = std::atan(out[i]) * computer_buffer[i];
                        }
                        derivator.compute(computer_buffer, out, nil);
                        for (int i = 0; i < data.size(); i++){
                            out[i] = std::abs(out[i]);
                        }
                        integrator.compute(out, computer_buffer, nil);
                        if constexpr (NP_DSP::CONFIG::debug){
                            computer_buffer.show(PlottingKind::Interpolate);
                        }
                        for (auto i = 0; i < data.size(); i++){
                            auto approx_answer_right = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer_right = approx_answer_right;
                            auto counter = 0;
                            while (approx_answer_right < std::numbers::pi){
                                counter++;
                                old_approx_answer_right = approx_answer_right;
                                approx_answer_right = computer_buffer.interpolate(i + counter) - computer_buffer.interpolate(i);
                            }
                            auto left_loss_right_edge = std::numbers::pi - old_approx_answer_right;
                            auto right_loss_right_edge = approx_answer_right - std::numbers::pi;
                            auto sum_loss_right_edge = left_loss_right_edge + right_loss_right_edge;
                            auto right_edge = static_cast<OutType::SampleType>(counter) - right_loss_right_edge/sum_loss_right_edge;
                            
                            auto approx_answer_left = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer_left = approx_answer_left;
                            counter = 0;
                            while (approx_answer_left < std::numbers::pi){
                                counter++;
                                old_approx_answer_left = approx_answer_left;
                                approx_answer_left = computer_buffer.interpolate(i)-computer_buffer.interpolate(i-counter);
                            }
                            auto left_loss_left_edge = std::numbers::pi - old_approx_answer_left;
                            auto right_loss_left_edge = approx_answer_left - std::numbers::pi;
                            auto sum_loss_left_edge = left_loss_left_edge + right_loss_left_edge;
                            auto left_edge = static_cast<OutType::SampleType>(counter) - right_loss_left_edge/sum_loss_left_edge; 

                            auto period = right_edge + left_edge;
                            out[i] = static_cast<OutType::SampleType>(1.0/period);
                        }
                    }
                }
            };

            export enum class ExtremumsBasedComputeInstFreqKind {Simple, Linear};

            export
            template<Signal DataT, Signal OutT, ExtremumsBasedComputeInstFreqKind compute_kind>
            struct ExtremumsBased{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = GENERAL::Nil;

                constexpr static ExtremumsBasedComputeInstFreqKind kind = compute_kind;
                constexpr static bool is_inst_freq_computer = true;

                void compute(const DataType & data, OutType & out, GENERAL::Nil nil){
                    std::vector<typename DataType::IdxType> extremums;
                    extremums.push_back(static_cast<typename DataType::IdxType>(0));
                    for (auto i = 1; i < data.size()-1; i++){
                        if ((data[i] >= data[i-1] &&
                            data[i] > data[i+1]) ||
                            (data[i] > data[i-1] &&
                            data[i] >= data[i+1]) ||
                            (data[i] <= data[i-1] &&
                             data[i] < data[i+1]) ||
                            (data[i] < data[i-1] &&
                             data[i] <= data[i+1]))
                        {
                            extremums.push_back(i);
                        }
                    }
                    extremums.push_back(static_cast<typename DataType::IdxType>(data.size()-1));
                    if(extremums.size() > 2){
                        extremums[0] = extremums[1] * 2 - extremums[2];
                        auto last = extremums.size() - 1;
                        extremums[last] = extremums[last - 1] * 2 - extremums[last - 2];
                        if (extremums[0] > 0){
                            extremums[0] = - extremums[1];
                        }
                        if (extremums[last] < extremums.size() - 1){
                            extremums[last] = data.size() + data.size() - extremums[last-1];
                        }
                    }
                    else{
                        auto last = extremums.size() - 1;
                        extremums[0] = - extremums[1];
                        extremums[last] = data.size() + data.size() - extremums[last-1];
                    }
                    
                    if constexpr (kind == ExtremumsBasedComputeInstFreqKind::Simple){
                        auto left = extremums[0];
                        auto right = extremums[1];
                        auto counter = 0;
                        for (auto i = 0; i < out.size(); i++){
                            if (i>right){
                                if (counter + 2 < extremums.size()){
                                    counter++;
                                }
                                left = extremums[counter];
                                right = extremums[counter+1];
                            }
                            out[i] = static_cast<OutType::SampleType>(0.5/(right - left));
                        }
                    }
                    else if constexpr (kind == ExtremumsBasedComputeInstFreqKind::Linear){
                        std::vector<std::pair<typename OutType::SampleType, typename OutType::SampleType>> points;
                        points.push_back({static_cast<OutType::SampleType>(0), static_cast<OutType::SampleType>(0.5/(extremums[1] - extremums[0]))});
                        for (auto i = 0; i < extremums.size()-1; i++){
                            points.push_back({static_cast<OutType::SampleType>((extremums[i+1] + extremums[i])/2.0),
                                              static_cast<OutType::SampleType>(0.5/(extremums[i+1] - extremums[i]))});
                        }
                        points.push_back({static_cast<OutType::SampleType>(data.size()-1), static_cast<OutType::SampleType>
                                    (0.5/(extremums[extremums.size()-1] - extremums[extremums.size()-2]))});

                        auto left = points[0];
                        auto right = points[1];
                        auto counter = 0;

                        for (auto i = 0; i < out.size(); i++){
                            if (i>right.first){
                                counter++;
                                left = points[counter];
                                right = points[counter+1];
                            }
                            out[i] = UTILITY_MATH::linearInterpolate(left, right,static_cast<OutType::SampleType>(i));
                        }
                    }
                    else{
                        std::unreachable();
                    }
                }
            };

            export
            template<Signal DataT, Signal OutT, /*SignalWrapper OptFn,*/
                    Integrator IntegratorT, Derivator DerivatorT,
            InstFreqDerivativeBasedKind kind>
            struct PeriodAndExtremumsBased {
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = GENERAL::Nil;

                using SampleType = DataType::SampleType;

                constexpr static bool is_inst_freq_computer = true;
                constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;

                float approx_order_coeff;

                //static_assert(OutType::is_writable == true);
                SampleType max_error = static_cast<typename DataType::SampleType>(1000000);

                PeriodAndExtremumsBased(IntegratorType integrator_in, DerivatorType derivator_in){
                    integrator = integrator_in;
                    derivator = derivator_in;
                }

                void compute(const DataType & data, OutType & out, AdditionalDataType & computer_buffer) {
                    GENERAL::Nil nil;
                    // compute extremums inst freq ->
                    // compute period based with external opt parameter equal const 1
                    // use fourier based approximator
                    // train it to approximate external opt parameter to minimize loss with extremums

                    std::vector<SampleType> extremums_freq_vec (data.size());
                    SimpleVecWrapper<SampleType> extremums_freq_wrapper(extremums_freq_vec);
                    GenericSignal<SimpleVecWrapper<SampleType>> extremums_freq(extremums_freq_wrapper);

                    ExtremumsBased<DataType, OutType, ExtremumsBasedComputeInstFreqKind::Linear> extremums_based;
                    extremums_based.compute(data, extremums_freq, nil);

                    std::vector<SampleType> external_opt_parametr_vector (data.size());
                    SimpleVecWrapper<SampleType> external_opt_parametr_wrapper(external_opt_parametr_vector);
                    GenericSignal<SimpleVecWrapper<SampleType>> external_opt_parametr(external_opt_parametr_wrapper);

                    DerivativeBasedWithExternalOptParametr<DataType, OutType, IntegratorType, DerivatorType, kind> 
                        inst_freq_computer(integrator, derivator);

                    for (auto i = 0; i < external_opt_parametr.size(); i++) {
                        external_opt_parametr[i] = static_cast<SampleType>(1.);
                    }
                    inst_freq_computer.compute(data, out, external_opt_parametr);

                    auto loss = [&](auto & approximator) {
                        approximator.is_actual = false;
                        std::vector<double> add_error = {};
                        for (auto i = 0; i < data.size(); i++){
                            external_opt_parametr[i] = approximator.compute(i);
                            if (external_opt_parametr[i] < 0.05){
                                add_error.push_back((0.05 - external_opt_parametr[i]));
                                external_opt_parametr[i] = 0.05;
                            }
                            else{
                                add_error.push_back(0.0);
                            }
                            
                        }
                        double accum = 0.0;
                        
                        //std::string mark = "show external_opt_parametr";
                        //IC(mark);
                        //external_opt_parametr.show(PlottingKind::Simple);
                        

                        inst_freq_computer.compute(data, out, external_opt_parametr);
                        
                        //mark = "show extr freq";
                        //IC(mark);
                        
                        //extremums_freq.show(PlottingKind::Simple);
                        

                        
                        //mark = "show out";
                        //IC(mark);
                        //out.show(PlottingKind::Simple);
                        
                        for(auto i = 0; i < extremums_freq.size(); i++) {

                            accum += add_error[i] + (extremums_freq[i] - out[i]) * (extremums_freq[i] - out[i]) / 1000; /// data[i];
                        }
                        IC(accum);
                        return accum;
                    };

                    auto stopPoint = [](auto losses_different, auto & approximator) {
                        if (losses_different > 0.00001){ //todo move precision to external parameter
                            return false;
                        }
                        else{
                            return true;
                        }
                    };
                    //APPROX::FourerSeriesBased
                    auto approximator = APPROX::FourierSeriesBasedPositive<DataType, decltype(loss), decltype(stopPoint)>
                        (loss, external_opt_parametr, stopPoint);
                        approximator.is_actual = false;
                    approximator.setpolynomsCount(data.size() / 2 * approx_order_coeff);
                    approximator.max_value = 10;
                    approximator.train();
                }
            };

            export
            template<Signal DataT, Signal OutT, /*SignalWrapper OptFn,*/
                    Integrator IntegratorT, Derivator DerivatorT,
            InstFreqDerivativeBasedKind kind>
            struct PeriodAndExtremumsBasedExternal {
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = OutT;

                using SampleType = DataType::SampleType;

                constexpr static bool is_inst_freq_computer = true;
                constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;

                float approx_order_coeff;

                //static_assert(OutType::is_writable == true);
                SampleType max_error = static_cast<typename DataType::SampleType>(1000000);

                PeriodAndExtremumsBasedExternal(IntegratorType integrator_in, DerivatorType derivator_in){
                    integrator = integrator_in;
                    derivator = derivator_in;
                }

                void compute(const DataType & data, OutType & out, AdditionalDataType & computer_buffer) {
                    GENERAL::Nil nil;
                    // compute extremums inst freq ->
                    // compute period based with external opt parameter equal const 1
                    // use fourier based approximator
                    // train it to approximate external opt parameter to minimize loss with extremums

                    std::vector<SampleType> extremums_freq_vec (data.size());
                    SimpleVecWrapper<SampleType> extremums_freq_wrapper(extremums_freq_vec);
                    GenericSignal<SimpleVecWrapper<SampleType>> extremums_freq(extremums_freq_wrapper);

                    ExtremumsBased<DataType, OutType, ExtremumsBasedComputeInstFreqKind::Linear> extremums_based;
                    extremums_based.compute(data, extremums_freq, nil);

                    std::vector<SampleType> external_opt_parametr_vector (data.size());
                    SimpleVecWrapper<SampleType> external_opt_parametr_wrapper(external_opt_parametr_vector);
                    GenericSignal<SimpleVecWrapper<SampleType>> external_opt_parametr(external_opt_parametr_wrapper);

                    DerivativeBased<DataType, OutType, IntegratorType, DerivatorType, kind> 
                        inst_freq_computer(integrator, derivator);

                    for (auto i = 0; i < external_opt_parametr.size(); i++) {
                        external_opt_parametr[i] = static_cast<SampleType>(1.);
                    }
                    inst_freq_computer.compute(data, out, computer_buffer);
                    
                    auto loss = [&](auto & approximator) {
                        approximator.is_actual = false;
                        for (auto i = 0; i < data.size(); i++){
                            external_opt_parametr[i] = approximator.compute(i);
                        }
                        double accum = 0.0;
                        
                        //std::string mark = "show external_opt_parametr";
                        //IC(mark);
                        //external_opt_parametr.show(PlottingKind::Simple);
                        
                        //std::vector<double> result = {};
                        for(auto i = 0; i < extremums_freq.size(); i++) {
                            if (out[i] != 0.0){
                                accum += (extremums_freq[i] - out[i] * external_opt_parametr[i])
                                    * (extremums_freq[i] - out[i] * external_opt_parametr[i]) / 1000 / out[i];
                            }
                            //result.push_back(out[i] * external_opt_parametr[i]);
                        }

                        //mark = "show result";
                        //IC(mark);
                        //matplot::plot(result);
                        //matplot::show();
                        IC(accum);
                        return accum;
                    };

                    auto stopPoint = [](auto losses_different, auto & approximator) {
                        if (losses_different > 0.000001){ //todo move precision to external parameter
                            return false;
                        }
                        else{
                            return true;
                        }
                    };

                    auto bySampleError = [&](auto & approximator, auto i){
                        return (extremums_freq[i] - out[i] * external_opt_parametr[i])
                                    * (extremums_freq[i] - out[i] * external_opt_parametr[i]) / 1000 / out[i];
                    };
                    //APPROX::FourerSeriesBased
                    auto approximator = APPROX::FourierSeriesBased<DataType, decltype(loss), decltype(stopPoint), APPROX::FSApproxKind::Simple, decltype(bySampleError)>
                        (loss, external_opt_parametr, stopPoint);
                    approximator.is_actual = false;
                    
                    approximator.setpolynomsCount(data.size() / 2 * approx_order_coeff);
                    approximator.max_value = 10;
                    approximator.train();
                    //approximator.fineTrainIter();
                    //approximator.train();
                    for (auto i = 0; i < data.size(); i++){
                        out[i] = out[i] * approximator.compute(i);
                    }
                    //todo
                }
            };
        }
    }
}