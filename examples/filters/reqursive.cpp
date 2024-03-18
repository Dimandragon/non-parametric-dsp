import signals;
import integrators;
import derivators;
import inst_freq_computers;
import phase_computers;
import inst_ampl_computers;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;
import filters;

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <ctime>

#include<icecream.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT inst_freq_buffer;
    SignalT inst_freq_mode_buffer;
    SignalT mode;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    //NP_DSP::ONE_D::
    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
            phase_computer;

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
            phase_computer_for_mode;

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer)> 
                    inst_freq_computer
                        (integrator, derivator, phase_computer);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer_for_mode)> 
                    inst_freq_computer_for_mode
                        (integrator, derivator, phase_computer_for_mode);
    
    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<double, decltype(integrator), decltype(derivator), decltype(inst_freq_computer),
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
            inst_ampl_computer(integrator, derivator, inst_freq_computer);

    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
        NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                non_opt_filter(integrator);

    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
        NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                non_opt_filter2(integrator);


    NP_DSP::ONE_D::FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
        decltype(non_opt_filter2)> cascade_filter(non_opt_filter, non_opt_filter2);

    NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilter<double, decltype(cascade_filter), decltype(inst_freq_computer),
        decltype(phase_computer), decltype(inst_freq_computer_for_mode), decltype(phase_computer_for_mode)> 
            opt_filter(cascade_filter, inst_freq_computer, phase_computer, 
                inst_freq_computer_for_mode, phase_computer_for_mode);

    NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChanges<double, decltype(integrator), decltype(derivator),
        decltype(cascade_filter), decltype(inst_freq_computer), decltype(phase_computer), decltype(inst_ampl_computer),
            decltype(inst_freq_computer_for_mode), decltype(phase_computer_for_mode)> filter
                (integrator, derivator, cascade_filter, inst_freq_computer, phase_computer, inst_ampl_computer,
                    inst_freq_computer_for_mode, phase_computer_for_mode);

    auto size = 50;
    for (auto i = 0; i < size; i++){
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
        inst_freq_buffer.base->vec->push_back(0);
        inst_freq_mode_buffer.base->vec->push_back(0);
        mode.base->vec->push_back(0);
    }
    SignalT aaaa;
    for (auto i = 0; i < 50; i++){
        aaaa.base->vec->push_back(std::rand());
    }

    inst_freq_computer.variability = 0.5;
    inst_freq_computer.compute(signal1, inst_freq_buffer, &compute_buffer);
    //opt_filter.compute(signal1, signal2, &inst_freq_buffer);
    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    /*std::string strLine = "This is a string!";
    std::cout << "warming up for tests ...\n";
    std::cout << strLine << std::endl;
    std::cout << "warmed up for tests ...\n";

    clock_t cpu_time_start;    // clock() time, CPU time used
    clock_t cpu_time_end;
    auto real_time_start = std::chrono::high_resolution_clock::now();   // real time, Wall clock time passed
    auto real_time_end = std::chrono::high_resolution_clock::now();

	//---------------------------------------------------------------------------------------------------
    
    cpu_time_start = clock();
    real_time_start = std::chrono::high_resolution_clock::now();*/
    filter.max_iter_number = 10000000000;
    aaaa.show(NP_DSP::ONE_D::PlottingKind::Simple);
    //for(int i = 0; i < 5; i++){
        filter.compute(signal1, signal2, &inst_freq_buffer);
        IC(filter.true_iter_number);
    //}

    /*std::cout << std::endl;

    cpu_time_end = clock();
    real_time_end = std::chrono::high_resolution_clock::now();
    std::cout << "iterator-based loop CPU time: " << 1000.0 * (cpu_time_end-cpu_time_start) / CLOCKS_PER_SEC << " ms\n"
         << "iterator-based loop real time: "<< std::chrono::duration<double, std::milli>(real_time_end-real_time_start).count() << " ms\n";
    */
}