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

namespace NP_DSP{
    namespace ONE_D{
        namespace FILTERS{
            export enum class InstFreqKind{Average, Double}; 
            export enum class FilteringType{DerivativeBased, ValueBased};

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

                NonOptPeriodBasedFilter(IntegratorT integrator){
                    //this->inst_freq = &inst_freq;
                    this->integrator = integrator;
                }

                void compute(const DataType & data, OutType & out, InstFreqType & inst_freq){
                    if constexpr (filtering_type_k == FilteringType::DerivativeBased){
                        if constexpr (inst_freq_kind == InstFreqKind::Average){
                            auto size_expr = [&](){
                                return data.size();
                            };

                            GENERAL::Nil nil;
                            auto val_expression = [&](DataType::IdxType idx){
                                return (data.interpolate(idx + 0.5/inst_freq.interpolate(idx)) - data.interpolate(idx - 0.5/inst_freq.interpolate(idx))) *
                                inst_freq.interpolate(idx);
                            };

                            ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType, 
                                decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                    expr_wrapper (val_expression, size_expr);

                            GenericSignal<decltype(expr_wrapper), false> expr_signal(expr_wrapper);
                            
                            NP_DSP::ONE_D::INTEGRATORS::Riman<decltype(expr_signal), OutType, 
                                NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator_new;
                            
                            //todo fix its big crucher
                            integrator_new.compute(expr_signal, out, nil);
                        }
                        else if constexpr (inst_freq_kind == InstFreqKind::Double){
                            /*auto val_expression = [=](DataType::IdxType idx){
                                return data.interpolate(idx + 0.5/inst_freq.interpolate(idx).forward) * inst_freq.interpolate(idx).forward
                                - data.interpolate(idx - 0.5/inst_freq.interpolate(idx).backward) * inst_freq.interpolate(idx).backward;
                            };
                            integrator.compute(ExpressionWrapper<typename DataType::SampleType,typename DataType::IdxType,
                                    decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                                       (val_expression, GENERAL::Nil{}, size_expr), out, {});*/
                            std::unreachable();
                        }
                        else{
                            std::unreachable();
                        }
                    }
                    else if constexpr (filtering_type_k == FilteringType::ValueBased){
                        if constexpr (inst_freq_kind == InstFreqKind::Average){
                            auto val_expression = [&](DataType::IdxType idx){
                                return (data.interpolate(idx + 0.25/inst_freq.interpolate(idx)) 
                                    + data.interpolate(idx - 0.25/inst_freq.interpolate(idx))) / 2;
                            };
                            for (auto i = 0; i < data.size(); i++){
                                out[i] = val_expression(i);
                            }
                        }
                        else if constexpr (inst_freq_kind == InstFreqKind::Double){
                            /*auto val_expression = [=](DataType::IdxType idx){
                                return data.interpolate(idx + 0.5/inst_freq.interpolate(idx).forward) * inst_freq.interpolate(idx).forward
                                - data.interpolate(idx - 0.5/inst_freq.interpolate(idx).backward) * inst_freq.interpolate(idx).backward;
                            };
                            integrator.compute(ExpressionWrapper<typename DataType::SampleType,typename DataType::IdxType,
                                    decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                                       (val_expression, GENERAL::Nil{}, size_expr), out, {});*/
                            std::unreachable();
                        }
                        else{
                            std::unreachable();
                        }
                    }
                    else{
                        std::unreachable();
                    }
                }
            };

            export 
            enum class InstFreqComputerKind {extremums_based, 
                arctg_extr_opt_momental, arctg_extr_opt_time_average,
                arctg_extr_opt_derive_average};

            export
            template<Signal DataT, Signal OutT, Signal InstFreqT, FilteringType filtering_type_k,
                    Integrator IntegratorT, Derivator DerivatorT, InstFreqComputerKind inst_freq_computer_k>
                    //INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind inst_freq_k>
            struct OptPeriodBasedFilter{
                using DataType = DataT;
                using OutType = OutT;
                using InstFreqType = InstFreqT;
                using AdditionalDataType = InstFreqType;
                IntegratorT integrator;
                DerivatorT derivator;
                NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> mode;
                NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> inst_freq_buffer2;

                double error_threshold = 0.0000001;
                using BuffT = decltype(mode);

                //constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
                constexpr static bool is_filter = true;
                NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<DataType, OutType, 
                        AdditionalDataType, filtering_type_k, 
                            IntegratorT, InstFreqKind::Average> * filter;
                
                auto * & inst_freq_computer(){
                    if constexpr (inst_freq_computer_k == InstFreqComputerKind::extremums_based){
                        static NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
                            <DataT, InstFreqT, 
                                NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear> 
                                    * inst_freq_computer;
                        return inst_freq_computer;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::arctg_extr_opt_momental){
                        static NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental> 
                                * inst_freq_computer;
                        return inst_freq_computer;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::arctg_extr_opt_time_average){
                        static NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage> 
                                * inst_freq_computer;
                        return inst_freq_computer;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::arctg_extr_opt_derive_average){
                        static NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage> 
                                * inst_freq_computer;
                        return inst_freq_computer;
                    }
                }

                auto * & inst_freq_computer_for_mode(){
                    if constexpr (inst_freq_computer_k == InstFreqComputerKind::extremums_based){
                        static NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
                            <BuffT, BuffT, 
                                NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear> 
                                    * inst_freq_computer_for_mode;
                        return inst_freq_computer_for_mode;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::arctg_extr_opt_momental){
                        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental> 
                                * inst_freq_computer_for_mode;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::arctg_extr_opt_time_average){
                        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage> 
                                * inst_freq_computer_for_mode;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::arctg_extr_opt_derive_average){
                        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage> 
                                * inst_freq_computer_for_mode;
                    }
                }
                /*
                
                
                
                */

                GENERAL::Nil nil;
                
                OptPeriodBasedFilter(IntegratorT integrator, DerivatorT derivator){
                    //this->inst_freq = &inst_freq;
                    this->integrator = integrator;
                    this->derivator = derivator;
                    filter = new NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter
                        <DataType, OutType, AdditionalDataType, filtering_type_k, 
                            IntegratorT, InstFreqKind::Average> (integrator);
                    
                    /*inst_freq_computer = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, 
                        OutT, IntegratorT, DerivatorT, inst_freq_k> (integrator, derivator); 

                    inst_freq_computer_for_mode = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT, 
                        DerivatorT, inst_freq_k> (integrator, derivator); */
                    if constexpr(inst_freq_computer_k == InstFreqComputerKind::extremums_based){
                        inst_freq_computer() = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
                            <DataT, InstFreqT, 
                                NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear> ;

                        inst_freq_computer_for_mode() = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
                            <BuffT, BuffT, 
                                NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear> ;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::arctg_extr_opt_momental){
                        inst_freq_computer() = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental> ;
                        
                        inst_freq_computer_for_mode() = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental> ;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::arctg_extr_opt_time_average){
                        inst_freq_computer() = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage> ;
                        
                        inst_freq_computer_for_mode() = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage> ;
                    }
                    else if constexpr (inst_freq_computer_k == InstFreqComputerKind::arctg_extr_opt_derive_average){
                        inst_freq_computer() = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage> ;
                        
                        inst_freq_computer_for_mode() = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT, 
                            DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage> ;
                    }
                }

                ~OptPeriodBasedFilter(){
                    delete filter;
                    delete inst_freq_computer();
                    delete inst_freq_computer_for_mode();
                }

                void computeIter(const DataType & data, OutType & out, InstFreqType & inst_freq_buffer){
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
                            mode[i] = (data[i] - out[i]);
                        }
                    }
                    if(inst_freq_buffer2.size() != inst_freq_buffer.size()){
                        inst_freq_buffer2.base->vec->clear();
                        for (int i = 0; i < inst_freq_buffer.size(); i++){
                            inst_freq_buffer2.base->vec->push_back(0);
                        }
                    }
                    out.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    mode.show(NP_DSP::ONE_D::PlottingKind::Simple);

                    //inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, out);
                    inst_freq_computer_for_mode()->compute(mode, inst_freq_buffer2, nil);

                    double error = UTILITY_MATH::signalsL2Distance
                        <double, InstFreqType, decltype(inst_freq_buffer2)>
                            (inst_freq_buffer, inst_freq_buffer2);
                    IC(error);
                    
                   

                    inst_freq_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    inst_freq_buffer2.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    if (error < error_threshold){
                        for (int i = 0; i < data.size(); i++){
                            out[i] = data[i] - mode[i];
                        }
                        return;
                    }
                    for (int i = 0; i < inst_freq_buffer.size(); i++){ //= i + std::rand() % (inst_freq_buffer.size()/5)){
                        inst_freq_buffer[i] = inst_freq_buffer[i] * (inst_freq_buffer2[i] + inst_freq_buffer[i] * 10.) / 11.0 / inst_freq_buffer[i];
                        //inst_freq_buffer.show();
                    }
                    computeIter(data, out, inst_freq_buffer);
                }

                void compute(const DataType & data, OutType & out, InstFreqType & inst_freq_buffer) {
                    //inst_freq_computer->compute(data, inst_freq_buffer, out);
                    inst_freq_computer()->compute(data, inst_freq_buffer, nil);
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
                            mode[i] = (data[i] - out[i]);
                        }
                    }
                    if(inst_freq_buffer2.size() != inst_freq_buffer.size()){
                        inst_freq_buffer2.base->vec->clear();
                        for (int i = 0; i < inst_freq_buffer.size(); i++){
                            inst_freq_buffer2.base->vec->push_back(0);
                        }
                    }
                    

                    //inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, out);
                    inst_freq_computer_for_mode()->compute(mode, inst_freq_buffer2, nil);

                    double error = UTILITY_MATH::signalsL2Distance
                        <double, InstFreqType, decltype(inst_freq_buffer2)>
                            (inst_freq_buffer, inst_freq_buffer2);
                    IC(error);
                    
                    mode.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    inst_freq_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    if (error < error_threshold){
                        for (int i = 0; i < data.size(); i++){
                            out[i] = data[i] - mode[i];
                        }
                        return;
                    }
                    for (int i = 0; i < inst_freq_buffer.size(); i++){
                        inst_freq_buffer[i] = inst_freq_buffer[i] * inst_freq_buffer2[i] / inst_freq_buffer[i];
                        //inst_freq_buffer.show();
                    }
                    computeIter(data, out, inst_freq_buffer);
                }
            };
/*
            export
            template<Signal DataT, Signal OutT, Signal InstFreqT, FilteringType filtering_type_k,
                    Integrator IntegratorT, Derivator DerivatorT, INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind inst_freq_k>
            struct OptFSPeriodBasedFilter{
                using DataType = DataT;
                using OutType = OutT;
                using InstFreqType = InstFreqT;
                using AdditionalDataType = InstFreqType;
                IntegratorT integrator;
                DerivatorT derivator;
                NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> mode;
                NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> inst_freq_buffer2;

                double error_threshold = 0.001;
                using BuffT = decltype(mode);

                //constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
                constexpr static bool is_filter = true;
                NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<DataType, OutType, 
                        AdditionalDataType, filtering_type_k, 
                            IntegratorT, InstFreqKind::Average> * filter;
                
                NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, OutT, IntegratorT, 
                    DerivatorT, inst_freq_k> * inst_freq_computer;
                
                NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT, 
                    DerivatorT, inst_freq_k> * inst_freq_computer_for_mode;

                OptPeriodBasedFilter(IntegratorT integrator, DerivatorT derivator){
                    //this->inst_freq = &inst_freq;
                    this->integrator = integrator;
                    this->derivator = derivator;
                    filter = new NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter
                        <DataType, OutType, AdditionalDataType, filtering_type_k, 
                            IntegratorT, InstFreqKind::Average> (integrator);
                    inst_freq_computer = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<DataT, InstFreqT, 
                        OutT, IntegratorT, DerivatorT, inst_freq_k> (integrator, derivator); 

                    inst_freq_computer_for_mode = new NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<BuffT, BuffT, OutT, IntegratorT, 
                        DerivatorT, inst_freq_k> (integrator, derivator); 
                }

                ~OptPeriodBasedFilter(){
                    delete filter;
                    delete inst_freq_computer;
                    delete inst_freq_computer_for_mode;
                }

                void computeIter(const DataType & data, OutType & out, InstFreqType & inst_freq_buffer){
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
                            mode[i] = (data[i] - out[i]);
                        }
                    }
                    if(inst_freq_buffer2.size() != inst_freq_buffer.size()){
                        inst_freq_buffer2.base->vec->clear();
                        for (int i = 0; i < inst_freq_buffer.size(); i++){
                            inst_freq_buffer2.base->vec->push_back(0);
                        }
                    }
                    out.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    mode.show(NP_DSP::ONE_D::PlottingKind::Simple);

                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, out);

                    double error = UTILITY_MATH::signalsL2Distance
                        <double, InstFreqType, decltype(inst_freq_buffer2)>
                            (inst_freq_buffer, inst_freq_buffer2);
                    IC(error);
                    
                   

                    inst_freq_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    if (error < error_threshold){
                        for (int i = 0; i < data.size(); i++){
                            out[i] = data[i] - mode[i];
                        }
                        return;
                    }
                    for (int i = 0; i < inst_freq_buffer.size(); i = i + std::rand() % 20){
                        inst_freq_buffer[i] = inst_freq_buffer[i] * (inst_freq_buffer2[i] + inst_freq_buffer[i] * 50.) / 51.0 / inst_freq_buffer[i];
                        //inst_freq_buffer.show();
                    }
                    computeIter(data, out, inst_freq_buffer);
                }

                void compute(const DataType & data, OutType & out, InstFreqType & inst_freq_buffer) {
                    inst_freq_computer->compute(data, inst_freq_buffer, out);
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
                            mode[i] = (data[i] - out[i]);
                        }
                    }
                    if(inst_freq_buffer2.size() != inst_freq_buffer.size()){
                        inst_freq_buffer2.base->vec->clear();
                        for (int i = 0; i < inst_freq_buffer.size(); i++){
                            inst_freq_buffer2.base->vec->push_back(0);
                        }
                    }
                    

                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, out);

                    double error = UTILITY_MATH::signalsL2Distance
                        <double, InstFreqType, decltype(inst_freq_buffer2)>
                            (inst_freq_buffer, inst_freq_buffer2);
                    IC(error);
                    
                    mode.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    inst_freq_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    if (error < error_threshold){
                        for (int i = 0; i < data.size(); i++){
                            out[i] = data[i] - mode[i];
                        }
                        return;
                    }
                    for (int i = 0; i < inst_freq_buffer.size(); i++){
                        inst_freq_buffer[i] = inst_freq_buffer[i] * inst_freq_buffer2[i] / inst_freq_buffer[i];
                        //inst_freq_buffer.show();
                    }
                    computeIter(data, out, inst_freq_buffer);
                }
            };*/
        }
        
    }
}