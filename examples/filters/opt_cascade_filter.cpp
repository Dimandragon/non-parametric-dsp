#include <signals.hpp>
#include <integrators.hpp>
#include <derivators.hpp>
#include <inst_freq_computers.hpp>
#include <phase_computers.hpp>
#include <inst_ampl_computers.hpp>
#include <cmath>
#include <vector>
#include <npdsp_concepts.hpp>
#include <string>
#include <filters.hpp>

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <ctime>

#include <icecream.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT inst_freq_buffer;
    SignalT mode;
    SignalT compute_buffer;
    

    

    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;

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

    NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
        NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                non_opt_filter(integrator);

    NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
        NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                non_opt_filter2(integrator);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
        NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage, decltype(phase_computer_for_mode)>
            inst_freq_computer_for_mode (integrator, derivator, phase_computer_for_mode);

    NP_DSP::ONE_D::FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
        decltype(non_opt_filter2)> cascade_filter(non_opt_filter, non_opt_filter2);

    NP_DSP::ONE_D::FILTERS::EXPERIMENTAL::OptPeriodBasedFilter<double, decltype(cascade_filter), decltype(inst_freq_computer),
        decltype(phase_computer), decltype(inst_freq_computer_for_mode), decltype(phase_computer_for_mode)> 
            filter(cascade_filter, inst_freq_computer, phase_computer, 
                inst_freq_computer_for_mode, phase_computer_for_mode);
    
    double counter_iters = 0.0;
    for(int size = 10; size < 2000; size++){
        signal1.base->vec->clear();
        signal2.base->vec->clear();
        compute_buffer.base->vec->clear();
        inst_freq_buffer.base->vec->clear();
        mode.base->vec->clear();
        for (int i = 0; i < size; i++){
            signal1.base->vec->push_back(std::rand());
            signal2.base->vec->push_back(0.);
            compute_buffer.base->vec->push_back(0.);
            inst_freq_buffer.base->vec->push_back(0.);
            mode.base->vec->push_back(0.);
        }
        inst_freq_computer.compute(signal1, inst_freq_buffer, &compute_buffer);

        filter.compute(signal1, signal2, &inst_freq_buffer);
        counter_iters += filter.true_iter_number;
        IC(filter.true_iter_number);
    }
    
    return 0;

    //cascade_filter.compute(signal1, signal2, &inst_freq_buffer);
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
    real_time_start = std::chrono::high_resolution_clock::now();
*/
    int size = 50;

        inst_freq_computer.compute(signal1, inst_freq_buffer, &compute_buffer);
        for(int j = 0; j < size; j++){
            signal2[j] = 0.0;
        }
        filter.compute(signal1, signal2, &inst_freq_buffer);
        IC(filter.true_iter_number);
/*
    std::cout << std::endl;
    cpu_time_end = clock();
    real_time_end = std::chrono::high_resolution_clock::now();
    std::cout << "iterator-based loop CPU time: " << 1000.0 * (cpu_time_end-cpu_time_start) / CLOCKS_PER_SEC << " ms\n"
         << "iterator-based loop real time: "<< std::chrono::duration<double, std::milli>(real_time_end-real_time_start).count() << " ms\n";
*/
    return 0;
    //IC(filter.error_old);
    
    for (int i = 0; i < size; i++){
        mode[i] = signal1[i] - signal2[i];
    }

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    mode.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;
}