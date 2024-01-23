module;

#include <icecream.hpp>

export module approximators;

import npdsp_concepts;
import npdsp_config;
import <vector>;
import <complex>;
import utility_math;
import <numbers>;
import npdsp_config;
import <string>;
import <matplot/matplot.h>;

namespace NP_DSP{
    namespace ONE_D{
        namespace APPROX{
            export
            template<Signal SignalT, typename LossFunc, typename StopPointFunc>
            struct FourierSeriesBased{
                using SignalType = SignalT;
                using SampleType = SignalT::SampleType;
                using Loss = LossFunc;
                using StopPoint = StopPointFunc;
                using IdxType = SignalT::IdxType;
                static constexpr bool is_signal_approximator = true;

                Loss * loss;
                StopPoint * stopPont;
                SignalType * signal;
                bool is_actual = false;


                std::vector<std::complex<SampleType>> fourier_series;
                std::vector<std::complex<SampleType>> approximated_data;

                SampleType max_value = 10000000000000000000000000.;


                FourierSeriesBased(Loss & lossFn, SignalType & signal_in, StopPointFunc & stop_point){
                    loss = &lossFn;
                    signal = &signal_in;
                    fourier_series = std::vector<std::complex<SampleType>>(signal_in.size());
                    approximated_data = std::vector<std::complex<SampleType>>(signal_in.size());
                    stopPont = &stop_point;
                }

                void computeData(){
                    is_actual = true;
                    UTILITY_MATH::ifftc2c(fourier_series, approximated_data);
                }

                std::complex<SampleType> computeSample(IdxType idx){
                    std::complex<SampleType> accum = {static_cast<SampleType>(0), static_cast<SampleType>(0)};
                    SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) / static_cast<SampleType>(fourier_series.size());
                    for (auto i = 0; i < fourier_series.size(); i++){
                        accum+= fourier_series[i]*std::complex<SampleType>{std::cos(w*i), std::sin(w*i)};
                    }
                    return accum;
                }

                SampleType compute(IdxType idx){
                    if(std::abs(static_cast<double>(idx) - static_cast<int64_t>(idx)) == 0){
                        if(is_actual){
                            return approximated_data[idx].real();
                        }
                        else{
                            computeData();
                            return approximated_data[idx].real();
                        }
                    }
                    else{
                        return computeSample(idx).real();
                    }
                }

                void train() {
                    for (auto i = 0; i < fourier_series.size(); i++) {
                        //auto trigonometric_sample = ONE_D::UTILITY_MATH::convertFSampleC2T(fourier_series[i]);
                        std::pair<SampleType, SampleType> trigonometric_sample;
                        trigonometric_sample.first = 10.;
                        //optimize
                        //optimize theta
                        if constexpr (CONFIG::debug){
                            std::string mark = "comopute fs cf";
                            IC(mark);
                            IC(i);
                        }
                        
                        if constexpr (CONFIG::debug){
                            std::string mark = "compute first losses";
                            IC(mark); 
                        }

                        SampleType check_loss[&](std::pair<T, T> data){
                            fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(data);
                            SampleType loss1 = (*loss)(*this);
                            fourier_series[i] = {fourier_series[i].imag2(), fourier_series[i].imag()};
                            SampleType loss2 = (*loss)(*this);
                            if (loss1<loss2)[
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(data);
                                return loss1;
                            ]
                            return loss2;
                        }

                        auto theta_min = static_cast<SampleType>(-std::numbers::pi/2.0 + 0.01);
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>({static_cast<SampleType>(10.), theta_min});
                        auto left_loss = (*loss)(*this);
                        auto theta_max = static_cast<SampleType>(std::numbers::pi/2.0 + 0.01);
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>({static_cast<SampleType>(10.), theta_max});
                        auto right_loss = (*loss)(*this);
                        auto theta_central = static_cast<SampleType>(0.0);
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>({static_cast<SampleType>(10.), theta_central});
                        auto central_loss = (*loss)(*this);

                        if constexpr (CONFIG::debug){
                            std::string mark = "end compute first losses";
                            IC(mark);
                        
                            IC(theta_min, left_loss, theta_central, central_loss, theta_max, right_loss);
                        }
                        
                        auto period_opt_iter = 0;

                        if constexpr (CONFIG::debug){
                            std::string mark = "start optimize period";
                            IC(mark); 
                        }
                        

                        while (!((*stopPont)(std::abs(right_loss - central_loss) + std::abs(left_loss - central_loss) , *this))) {
                            if constexpr (CONFIG::debug){
                                IC(period_opt_iter,left_loss, central_loss, right_loss, theta_min, theta_central, theta_max);
                                IC(right_loss-central_loss, left_loss-central_loss, right_loss-left_loss);
                                IC(theta_central - theta_min, theta_max-theta_central);
                            }
                            
                            auto left_diff = left_loss - central_loss;
                            auto right_diff = right_loss - central_loss;
                            if (left_diff > 0 && right_diff <=0) {
                                left_loss = central_loss;
                                theta_min = theta_central;
                                theta_central = (theta_min + theta_max) / 2;
                                trigonometric_sample.second = theta_central;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                central_loss = (*loss)(*this);
                            }
                            else if (left_diff <= 0 && right_diff > 0) {
                                right_loss = central_loss;
                                theta_max = theta_central;
                                theta_central = (theta_min + theta_max) / 2;
                                trigonometric_sample.second = theta_central;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                central_loss = (*loss)(*this);
                            }
                            else if (left_diff > right_diff) {
                                //left branch is higher
                                auto theta_left_avg = (theta_min + theta_central) / 2;
                                trigonometric_sample.second = theta_left_avg;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                auto new_left_loss = (*loss)(*this);
                                if(new_left_loss < left_loss && new_left_loss > central_loss) {
                                    left_loss = new_left_loss;
                                    theta_min = theta_left_avg;
                                    theta_central = (theta_min + theta_max) / 2;
                                    trigonometric_sample.second = theta_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                    central_loss = (*loss)(*this);
                                }
                                else {
                                    right_loss = central_loss;
                                    theta_max = theta_central;
                                    theta_central = theta_left_avg;
                                    central_loss = new_left_loss;
                                }
                            }
                            else if (left_diff <= right_diff){
                                //right branch is higher
                                auto theta_right_avg = (theta_max + theta_central) / 2;
                                trigonometric_sample.second = theta_right_avg;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                auto new_right_loss = (*loss)(*this);
                                if (new_right_loss < right_loss && new_right_loss > central_loss) {
                                    right_loss = new_right_loss;
                                    theta_max = theta_right_avg;
                                    theta_central = (theta_max + theta_min) / 2;
                                    trigonometric_sample.second = theta_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                    central_loss = (*loss)(*this);
                                }
                                else {
                                    left_loss = central_loss;
                                    theta_min = theta_central;
                                    central_loss = new_right_loss;
                                    theta_central = theta_right_avg;
                                }
                            }
                            period_opt_iter++;
                        }
                        trigonometric_sample.second = theta_central;
                        

                        if constexpr (CONFIG::debug){
                            std::string mark = "start optimize ampl";
                            IC(mark);   
                        }
                        
                        SampleType ampl_left = 0.0;
                        SampleType ampl_right = max_value;
                        SampleType ampl_central = (ampl_left + ampl_right) / 2;
                        trigonometric_sample.first = ampl_left;
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                        auto loss_left = (*loss)(*this);
                        trigonometric_sample.first = ampl_right;
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                        auto loss_right = (*loss)(*this);
                        trigonometric_sample.first = ampl_central;
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                        auto loss_central = (*loss)(*this);
                        auto ampl_opt_iter = 0;
                        //while(!(*stopPont)(std::abs(ampl_right - ampl_left), *this)) {
                        auto errors_counter = 0; 
                        while (!((*stopPont)(std::abs(loss_right - loss_central) + std::abs(loss_left - loss_central), *this))) {
                            if constexpr (CONFIG::debug){
                                IC(ampl_opt_iter, loss_left, loss_central, loss_right, ampl_left, ampl_central, ampl_right);
                                IC(loss_right-loss_central, loss_left-loss_central, loss_right-loss_left);
                                IC(ampl_central - ampl_left, ampl_right - ampl_central);
                            }
                            
                            auto left_diff = loss_left-loss_central;
                            auto right_diff = loss_right-loss_central;
                            if(left_diff > 0 && right_diff <=0) {
                                loss_left = loss_central;
                                ampl_left = ampl_central;
                                ampl_central = (ampl_left + ampl_right) / 2;
                                trigonometric_sample.first = ampl_central;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                loss_central = (*loss)(*this);
                            }
                            else if (left_diff <= 0 && right_diff > 0) {
                                loss_right = loss_central;
                                ampl_right = ampl_central;
                                ampl_central = (ampl_left + ampl_right) / 2;
                                trigonometric_sample.first = ampl_central;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                loss_central = (*loss)(*this);
                            }
                            else if (left_diff > right_diff) {
                                //left branch is higher
                                auto ampl_left_avg = (ampl_left + ampl_central) / 2;
                                trigonometric_sample.first = ampl_left_avg;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                auto new_loss_left = (*loss)(*this);
                                if (new_loss_left < loss_left && new_loss_left > central_loss) {
                                    loss_left = new_loss_left;
                                    ampl_left = ampl_left_avg;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                    loss_central = (*loss)(*this);
                                }
                                else {
                                    loss_right = loss_central;
                                    ampl_right = ampl_central;
                                    ampl_central = ampl_left_avg;
                                    loss_central = new_loss_left;
                                }
                            }
                            else if (left_diff <= right_diff) {
                                auto ampl_right_avg = (ampl_right + ampl_central) / 2;
                                trigonometric_sample.first = ampl_right_avg;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                auto new_loss_right = (*loss)(*this);
                                if (new_loss_right < loss_right && new_loss_right > loss_central) {
                                    loss_right = new_loss_right;
                                    ampl_right = ampl_right_avg;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                                    loss_central = (*loss)(*this);
                                }
                                else {
                                    loss_left = loss_central;
                                    ampl_left = ampl_central;
                                    loss_central = new_loss_right;
                                    ampl_central = ampl_right_avg;
                                }
                            }
                            //auto left_diff = loss_left-loss_central;
                            //auto right_diff = loss_right-loss_central;
                            if (left_diff == (loss_left-loss_central)){
                                if (right_diff == (loss_right-loss_central)){
                                    //std::string point_mark = "aaaa fuck x3";
                                    //IC(point_mark);
                                    if constexpr (CONFIG::debug){
                                        IC(left_diff, loss_left-loss_central, right_diff, loss_right-loss_central);
                                    }
                                    
                                    errors_counter++;
                                    if constexpr (CONFIG::debug){
                                        IC(errors_counter);
                                    }
                                    
                                    if (errors_counter > 10){
                                        if constexpr (CONFIG::debug){
                                            IC("aaaa fuck x4");
                                            for(;;){}
                                        }
                                        
                                    }
                                }
                                else{
                                    //std::string point_mark = "aaaa fuck x2";
                                    //IC(point_mark);
                                    if constexpr (CONFIG::debug){
                                        IC(right_diff == (loss_right-loss_central));
                                        IC(left_diff, loss_left-loss_central, right_diff, loss_right-loss_central);
                                    }
                                }
                            }
                            else{
                                if constexpr (CONFIG::debug){
                                    std::string point_mark = "aaaa fuck";
                                    IC(point_mark);
                                    IC(left_diff, loss_left-loss_central, right_diff, loss_right-loss_central);
                                }
                            }
                            ampl_opt_iter++;
                        }
                        trigonometric_sample.first = ampl_central;
                        
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C<SampleType>(trigonometric_sample);
                        show(PlottingKind::Simple);
                    }
                }


                void show(PlottingKind kind){
                    is_actual = false;
                    if (kind == PlottingKind::Simple){
                        std::vector<SampleType> plotting_data = {};
                        for (auto i = 0; i < signal->size(); i++){
                            plotting_data.push_back(compute(i));
                        }
                        matplot::plot(plotting_data);
                        matplot::show();
                    }
                }

                void show(PlottingKind kind, const std::string & filename, const std::string & format){
                    is_actual = false;
                    if (kind == PlottingKind::Simple){
                        std::vector<SampleType> plotting_data = {};
                        for (auto i = 0; i < signal->size(); i++){
                            plotting_data.push_back(compute(i));
                        }
                        matplot::plot(plotting_data);
                        matplot::show();
                        matplot::save(filename, format);
                    }
                }

                void show(PlottingKind kind, const std::string & filename){
                    is_actual = false;
                    if (kind == PlottingKind::Simple){
                        std::vector<SampleType> plotting_data = {};
                        for (auto i = 0; i < signal->size(); i++){
                            plotting_data.push_back(compute(i));
                        }
                        matplot::plot(plotting_data);
                        matplot::show();
                        matplot::save(filename);
                    }
                }
            };
        }
    }
}