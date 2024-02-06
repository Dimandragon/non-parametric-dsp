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
            export enum class FSApproxKind {Simple, Positive};
            export
            template<Signal SignalT, typename LossFunc, typename StopPointFunc, FSApproxKind kind_v, typename BySampleLoss>
            struct FourierSeriesBased{
                using SignalType = SignalT;
                using SampleType = SignalT::SampleType;
                using Loss = LossFunc;
                using StopPoint = StopPointFunc;
                using IdxType = SignalT::IdxType;
                static constexpr bool is_signal_approximator = true;

                static constexpr FSApproxKind kind = kind_v;

                Loss * loss;
                BySampleLoss * bySampleLoss = NULL;

                StopPoint * stopPont;
                SignalType * signal;
                bool is_actual = false;
                int polynoms_count;

                int polynoms_count_on_tile;
                int tile_size;

                std::vector<std::complex<SampleType>> fourier_series;
                std::vector<std::complex<SampleType>> approximated_data;

                SampleType max_value = 10000000000000000000000000.;


                FourierSeriesBased(Loss & lossFn, SignalType & signal_in, StopPointFunc & stop_point){
                    loss = &lossFn;
                    signal = &signal_in;
                    fourier_series = std::vector<std::complex<SampleType>>(signal_in.size());
                    approximated_data = std::vector<std::complex<SampleType>>(signal_in.size());
                    stopPont = &stop_point;
                    polynoms_count = signal_in.size() / 2;
                    polynoms_count_on_tile = signal_in.size() / 2;
                    tile_size = signal_in.size();
                }

                size_t mirrorIdx(size_t i){
                    return i / tile_size * tile_size + i / tile_size * tile_size + tile_size - i;
                }

                void applyMirror(size_t i){
                    if (i % tile_size != 0){
                        auto mirror_idx = mirrorIdx(i);
                        fourier_series[mirror_idx] = {fourier_series[i].real() / 2.0, -fourier_series[i].imag() / 2.0};
                        fourier_series[i] = fourier_series[i] / 2.0;
                    }
                }

                void computeData(){
                    is_actual = true;
                    for (auto i = 0; i < approximated_data.size() / tile_size; i++){
                        auto pad = i * tile_size;
                        UTILITY_MATH::ifftc2c(fourier_series, approximated_data, tile_size, pad);
                    }
                    auto pad = approximated_data.size() / tile_size * tile_size;
                    UTILITY_MATH::ifftc2c(fourier_series, approximated_data, approximated_data.size() - pad, pad);
                }

                void computeTile(size_t idx){
                    auto i = idx/tile_size;
                    auto pad = i * tile_size;
                    if (approximated_data.size() > pad + tile_size){
                        UTILITY_MATH::ifftc2c(fourier_series, approximated_data, tile_size, pad);
                    }
                    else{
                        UTILITY_MATH::ifftc2c(fourier_series, approximated_data, approximated_data.size() - pad, pad);
                    }
                }

                void setpolynomsCount(int count){
                    polynoms_count_on_tile = count;
                }

                std::complex<SampleType> computeSample(IdxType idx){
                    std::complex<SampleType> accum = {static_cast<SampleType>(0), static_cast<SampleType>(0)};
                    
                    if (idx<=approximated_data.size() && idx >= 0){
                        auto tile_first_idx = idx / tile_size * tile_size;
                        if (tile_first_idx + tile_size < approximated_data.size()){
                            SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) / static_cast<SampleType>(tile_size);
                            for (auto i = 0; i < tile_size; i++){
                                accum+= fourier_series[i + tile_first_idx]*std::complex<SampleType>{std::cos(w*i), std::sin(w*i)};
                            }
                        }
                        else{
                            SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) / 
                                static_cast<SampleType>(approximated_data.size() - tile_first_idx);
                            for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++){
                                accum+= fourier_series[i + tile_first_idx]*std::complex<SampleType>{std::cos(w*i), std::sin(w*i)};
                            }
                        }
                    }
                    else if (idx > approximated_data.size()){
                        auto tile_first_idx = approximated_data.size() / tile_size * tile_size;
                        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) / 
                                static_cast<SampleType>(approximated_data.size() - tile_first_idx);
                        for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++){
                            accum+= fourier_series[i + tile_first_idx]*std::complex<SampleType>{std::cos(w*i), std::sin(w*i)};
                        }
                    }
                    else{
                        if (approximated_data.size() > tile_size){
                            SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) / 
                                static_cast<SampleType>(tile_size);
                            for (auto i = 0; i < tile_size; i++){
                                accum+= fourier_series[i]*std::complex<SampleType>{std::cos(w*i), std::sin(w*i)};
                            }
                        }
                        else{
                            SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) / 
                                static_cast<SampleType>(approximated_data.size());
                            for (auto i = 0; i < approximated_data.size(); i++){
                                accum+= fourier_series[i]*std::complex<SampleType>{std::cos(w*i), std::sin(w*i)};
                            }
                        }
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

                std::complex<SampleType> computeComplex(IdxType idx){
                    if(std::abs(static_cast<double>(idx) - static_cast<int64_t>(idx)) == 0){
                        if(is_actual){
                            return approximated_data[idx];
                        }
                        else{
                            computeData();
                            return approximated_data[idx];
                        }
                    }
                    else{
                        return computeSample(idx);
                    }
                }

                void fineTrainIter(){
                    computeData();
                    //todo for tiles and for positive
                    for (auto i = 0; i < fourier_series.size(); i++){
                        if (i % tile_size > tile_size/2)
                        {
                            //auto idx_mirror = mirrorIdx(i);
                            //fourier_series[i] = {fourier_series[idx_mirror].real() / 2.0, -fourier_series[idx_mirror].imag() / 2.0};
                            //fourier_series[idx_mirror] = fourier_series[idx_mirror] / 2.0;
                            continue;
                        }
                        if (i % tile_size >= polynoms_count_on_tile){
                            continue;
                        }


                        auto check_loss = [&](std::pair<SampleType, SampleType> data){
                            auto complex_sample = UTILITY_MATH::convertFSampleT2C<SampleType>(data);
                            fourier_series[i] = complex_sample;
                            applyMirror(i);
                            computeTile(i);
                            SampleType loss1 = 0.;
                            if (bySampleLoss) {
                                size_t pad = i / tile_size * tile_size;
                                for (auto idx = 0; idx < tile_size; idx++){
                                    if (idx + pad >= approximated_data.size()){
                                        break;
                                    }
                                    loss1 += (*bySampleLoss)(*this, idx + pad);
                                }
                            }
                            else{
                                loss1 = (*loss)(*this);
                            }
                            fourier_series[i] = {complex_sample.real(), -complex_sample.imag()};
                            applyMirror(i);
                            computeTile(i);
                            SampleType loss2 = 0.;
                            if (bySampleLoss) {
                                size_t pad = i / tile_size * tile_size;
                                for (auto idx = 0; idx < tile_size; idx++){
                                    if (idx + pad >= approximated_data.size()){
                                        break;
                                    }
                                    loss2 += (*bySampleLoss)(*this, idx + pad);
                                }
                            }
                            else{
                                loss2 = (*loss)(*this);
                            }
                            fourier_series[i] = {-complex_sample.real(), complex_sample.imag()};
                            applyMirror(i);
                            computeTile(i);
                            SampleType loss3 = 0.;
                            if (bySampleLoss) {
                                size_t pad = i / tile_size * tile_size;
                                for (auto idx = 0; idx < tile_size; idx++){
                                    if (idx + pad >= approximated_data.size()){
                                        break;
                                    }
                                    loss3 += (*bySampleLoss)(*this, idx + pad);
                                }
                            }
                            else{
                                loss3 = (*loss)(*this);
                            }
                            fourier_series[i] = {-complex_sample.real(), -complex_sample.imag()};
                            applyMirror(i);
                            computeTile(i);
                            SampleType loss4 = 0.;
                            if (bySampleLoss) {
                                size_t pad = i / tile_size * tile_size;
                                for (auto idx = 0; idx < tile_size; idx++){
                                    if (idx + pad >= approximated_data.size()){
                                        break;
                                    }
                                    loss4 += (*bySampleLoss)(*this, idx + pad);
                                }
                            }
                            else{
                                loss4 = (*loss)(*this);
                            }

                            if (loss1 <= loss2 && loss1 <= loss3 && loss1 <= loss4){
                                fourier_series[i] = complex_sample;
                                applyMirror(i);
                                computeTile(i);
                                return loss1;
                            }
                            if (loss2 <= loss1 && loss2 <= loss3 && loss2 <= loss4){
                                fourier_series[i] = {complex_sample.real(), -complex_sample.imag()};
                                applyMirror(i);
                                computeTile(i);
                                return loss2;
                            }
                            if (loss3 <= loss2 && loss3 <= loss1 && loss3 <= loss4){
                                fourier_series[i] = {-complex_sample.real(), complex_sample.imag()};
                                applyMirror(i);
                                computeTile(i);
                                return loss3;
                            }

                            fourier_series[i] = {-complex_sample.real(), -complex_sample.imag()};
                            applyMirror(i);
                            computeTile(i);
                            return loss4;
                        };

                        auto trigonometric_sample = ONE_D::UTILITY_MATH::convertFSampleC2T(fourier_series[i]);

                        //auto theta_right = trigonometric_sample.second * 1.3;
                        auto theta_max = trigonometric_sample.second + std::numbers::pi/10.0;
                        auto theta_central = trigonometric_sample.second;
                        //auto theta_left = trigonometric_sample.second * 0.7;
                        auto theta_min = trigonometric_sample.second - std::numbers::pi/10.0;
                        if (theta_max > std::numbers::pi/2.0){
                            theta_max = std::numbers::pi/2.0;
                        }
                        if (theta_min < -std::numbers::pi/2.0){
                            theta_min = -std::numbers::pi/2.0;
                        }

                        auto ampl = trigonometric_sample.first;

                        auto left_loss = check_loss({static_cast<SampleType>(1.), theta_min});
                        auto central_loss = check_loss({static_cast<SampleType>(1.), theta_central});
                        auto right_loss = check_loss({static_cast<SampleType>(1.), theta_max});

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
                                central_loss = check_loss(trigonometric_sample);
                            }
                            else if (left_diff <= 0 && right_diff > 0) {
                                right_loss = central_loss;
                                theta_max = theta_central;
                                theta_central = (theta_min + theta_max) / 2;
                                trigonometric_sample.second = theta_central;
                                central_loss = check_loss(trigonometric_sample);
                            }
                            else if (left_diff > right_diff) {
                                //left branch is higher
                                auto theta_left_avg = (theta_min + theta_central) / 2;
                                trigonometric_sample.second = theta_left_avg;
                                auto new_left_loss = check_loss(trigonometric_sample);
                                if(new_left_loss < left_loss && new_left_loss > central_loss) {
                                    left_loss = new_left_loss;
                                    theta_min = theta_left_avg;
                                    theta_central = (theta_min + theta_max) / 2;
                                    trigonometric_sample.second = theta_central;
                                    central_loss = check_loss(trigonometric_sample);
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
                                auto new_right_loss = check_loss(trigonometric_sample);
                                if (new_right_loss < right_loss && new_right_loss > central_loss) {
                                    right_loss = new_right_loss;
                                    theta_max = theta_right_avg;
                                    theta_central = (theta_max + theta_min) / 2;
                                    trigonometric_sample.second = theta_central;
                                    central_loss = check_loss(trigonometric_sample);
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
                        SampleType ampl_twenty_procent = ampl * 0.2;
                        SampleType ampl_left = ampl;
                        SampleType ampl_right = ampl;
                        SampleType ampl_central = ampl;
                        if (ampl_twenty_procent < 2.){
                            ampl_left-= 2.;
                            ampl_right+=2.;
                        }
                        else{
                            ampl_left-=ampl_twenty_procent;
                            ampl_right+=ampl_twenty_procent;
                        }
                        if (ampl_left < 0.){
                            ampl_left = 0.;
                        }
                        trigonometric_sample.first = ampl_left;
                        auto loss_left = check_loss(trigonometric_sample);
                        trigonometric_sample.first = ampl_right;
                        auto loss_right = check_loss(trigonometric_sample);
                        trigonometric_sample.first = ampl_central;
                        auto loss_central = check_loss(trigonometric_sample);
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
                                std::string point_mark = "ampl 1 branch";
                                IC(point_mark);
                                loss_left = loss_central;
                                ampl_left = ampl_central;
                                ampl_central = (ampl_left + ampl_right) / 2;
                                trigonometric_sample.first = ampl_central;
                                loss_central = check_loss(trigonometric_sample);
                            }
                            else if (left_diff <= 0 && right_diff > 0) {
                                std::string point_mark = "ampl 2 branch";
                                IC(point_mark);
                                loss_right = loss_central;
                                ampl_right = ampl_central;
                                ampl_central = (ampl_left + ampl_right) / 2;
                                trigonometric_sample.first = ampl_central;
                                loss_central = check_loss(trigonometric_sample);
                            }
                            else if (left_diff > right_diff) {
                                std::string point_mark = "ampl 3 branch";
                                IC(point_mark);
                                //left branch is higher
                                auto ampl_left_avg = (ampl_left + ampl_central) / 2;
                                trigonometric_sample.first = ampl_left_avg;
                                auto new_loss_left = check_loss(trigonometric_sample);
                                if (new_loss_left < loss_left && new_loss_left > central_loss) {
                                    std::string point_mark = "ampl 3.1 branch";
                                    IC(point_mark);
                                    loss_left = new_loss_left;
                                    ampl_left = ampl_left_avg;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    loss_central = check_loss(trigonometric_sample);
                                }
                                else {
                                    std::string point_mark = "ampl 3.2 branch";
                                    IC(point_mark);
                                    ampl_right = (ampl_central+ampl_right) / 2;
                                    trigonometric_sample.first = ampl_right;
                                    loss_right = check_loss(trigonometric_sample);
                                    ampl_central = (ampl_right + ampl_left) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    loss_central = check_loss(trigonometric_sample);
                                }
                            }
                            else if (left_diff <= right_diff) {
                                std::string point_mark = "ampl 4 branch";
                                IC(point_mark);
                                auto ampl_right_avg = (ampl_right + ampl_central) / 2;
                                trigonometric_sample.first = ampl_right_avg;
                                auto new_loss_right = check_loss(trigonometric_sample);
                                if (new_loss_right < loss_right && new_loss_right > loss_central) {
                                    std::string point_mark = "ampl 4.1 branch";
                                    IC(point_mark);
                                    loss_right = new_loss_right;
                                    ampl_right = ampl_right_avg;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    loss_central = check_loss(trigonometric_sample);
                                }
                                else {
                                    std::string point_mark = "ampl 4.2 branch";
                                    IC(point_mark, i);
                                    ampl_left = (ampl_central + ampl_left)/2;
                                    trigonometric_sample.first = ampl_left;
                                    //IC(point_mark);
                                    loss_left = check_loss(trigonometric_sample);
                                    ampl_central = (ampl_right + ampl_left) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    loss_central = check_loss(trigonometric_sample);
                                    //IC(ampl_left, ampl_central, ampl_right, loss_left, loss_central, loss_right);
                                    //IC(point_mark);
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
                                            //IC("aaaa fuck x4");
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
                            if (ampl_central - ampl_left < 0.000001 || ampl_right - ampl_central < 0.000001){
                                break;
                            }
                        }
                        if (loss_central <= loss_right && loss_central <= loss_left){
                            trigonometric_sample.first = ampl_central;
                        }
                        else if (loss_right<loss_left){
                            trigonometric_sample.first = ampl_right;
                        }
                        else{
                            trigonometric_sample.first = ampl_left;
                        }
                        
                        
                        check_loss(trigonometric_sample);
                        //show(PlottingKind::Simple);
                    }
                }

                void train() {
                    computeData();
                    for (auto i = 0; i < fourier_series.size(); i++) {
                        // сделать это по тайлам
                        if (i % tile_size > tile_size/2)
                        {
                            //auto idx_mirror = i / tile_size * tile_size + i / tile_size * tile_size + tile_size - i;
                            //fourier_series[i] = {fourier_series[idx_mirror].real() / 2.0, -fourier_series[idx_mirror].imag() / 2.0};
                            //fourier_series[idx_mirror] = fourier_series[idx_mirror] / 2.0;
                            continue;
                        }
                        if (i % tile_size >= polynoms_count_on_tile){
                            continue;
                        }
                        std::pair<SampleType, SampleType> trigonometric_sample;
                        trigonometric_sample.first = 1.;
                        if constexpr (CONFIG::debug){
                            std::string mark = "comopute fs cf";
                            IC(mark);
                            IC(i);
                        }
                        
                        if constexpr (CONFIG::debug){
                            std::string mark = "compute first losses";
                            IC(mark); 
                        }

                        auto check_loss = [&](std::pair<SampleType, SampleType> data){
                            auto complex_sample = UTILITY_MATH::convertFSampleT2C<SampleType>(data);
                            fourier_series[i] = complex_sample;
                            SampleType loss1 = (*loss)(*this);
                            fourier_series[i] = {complex_sample.real(), -complex_sample.imag()};
                            SampleType loss2 = (*loss)(*this);
                            fourier_series[i] = {-complex_sample.real(), complex_sample.imag()};
                            SampleType loss3 = (*loss)(*this);
                            fourier_series[i] = {-complex_sample.real(), -complex_sample.imag()};
                            SampleType loss4 = (*loss)(*this);

                            if (loss1 <= loss2 && loss1 <= loss3 && loss1 <= loss4){
                                fourier_series[i] = complex_sample;
                                return loss1;
                            }
                            else if (loss2 <= loss1 && loss2 <= loss3 && loss2 <= loss4){
                                fourier_series[i] = {complex_sample.real(), -complex_sample.imag()};
                                return loss2;
                            }
                            else if (loss3 <= loss2 && loss3 <= loss1 && loss3 <= loss4){
                                fourier_series[i] = {-complex_sample.real(), complex_sample.imag()};
                                return loss3;
                            }
                            else{
                                fourier_series[i] = {-complex_sample.real(), -complex_sample.imag()};
                                return loss4;
                            }
                        };

                        auto theta_min = static_cast<SampleType>(-std::numbers::pi/2.0 + 0.01);
                        auto theta_max = static_cast<SampleType>(std::numbers::pi/2.0 + 0.01);
                        auto theta_central = static_cast<SampleType>(0.0);
                        auto max_ampl = 1.;
                        if constexpr (kind == FSApproxKind::Positive){
                            if (i % tile_size != 0){
                                auto thr_sample = UTILITY_MATH::convertFSampleC2T<SampleType>(fourier_series[0]);
                            
                                max_ampl = thr_sample.first * std::cos(thr_sample.second);
                                //IC(max_ampl);
                            }
                        }
                        auto left_loss = check_loss({static_cast<SampleType>(max_ampl), theta_min});
                        auto right_loss = check_loss({static_cast<SampleType>(max_ampl), theta_max});
                        auto central_loss = check_loss({static_cast<SampleType>(max_ampl), theta_central});

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
                                central_loss = check_loss(trigonometric_sample);
                            }
                            else if (left_diff <= 0 && right_diff > 0) {
                                right_loss = central_loss;
                                theta_max = theta_central;
                                theta_central = (theta_min + theta_max) / 2;
                                trigonometric_sample.second = theta_central;
                                central_loss = check_loss(trigonometric_sample);
                            }
                            else if (left_diff > right_diff) {
                                //left branch is higher
                                auto theta_left_avg = (theta_min + theta_central) / 2;
                                trigonometric_sample.second = theta_left_avg;
                                auto new_left_loss = check_loss(trigonometric_sample);
                                if(new_left_loss < left_loss && new_left_loss > central_loss) {
                                    left_loss = new_left_loss;
                                    theta_min = theta_left_avg;
                                    theta_central = (theta_min + theta_max) / 2;
                                    trigonometric_sample.second = theta_central;
                                    central_loss = check_loss(trigonometric_sample);
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
                                auto new_right_loss = check_loss(trigonometric_sample);
                                if (new_right_loss < right_loss && new_right_loss > central_loss) {
                                    right_loss = new_right_loss;
                                    theta_max = theta_right_avg;
                                    theta_central = (theta_max + theta_min) / 2;
                                    trigonometric_sample.second = theta_central;
                                    central_loss = check_loss(trigonometric_sample);
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
                        if constexpr (kind == FSApproxKind::Positive){
                            if (i!=0){
                                ampl_right = max_ampl;
                            }
                        }
                        SampleType ampl_central = (ampl_left + ampl_right) / 2;
                        trigonometric_sample.first = ampl_left;
                        auto loss_left = check_loss(trigonometric_sample);
                        trigonometric_sample.first = ampl_right;
                        auto loss_right = check_loss(trigonometric_sample);
                        trigonometric_sample.first = ampl_central;
                        auto loss_central = check_loss(trigonometric_sample);
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
                                std::string point_mark = "ampl 1 branch";
                                IC(point_mark);
                                loss_left = loss_central;
                                ampl_left = ampl_central;
                                ampl_central = (ampl_left + ampl_right) / 2;
                                trigonometric_sample.first = ampl_central;
                                loss_central = check_loss(trigonometric_sample);
                            }
                            else if (left_diff <= 0 && right_diff > 0) {
                                std::string point_mark = "ampl 2 branch";
                                IC(point_mark);
                                loss_right = loss_central;
                                ampl_right = ampl_central;
                                ampl_central = (ampl_left + ampl_right) / 2;
                                trigonometric_sample.first = ampl_central;
                                loss_central = check_loss(trigonometric_sample);
                            }
                            else if (left_diff > right_diff) {
                                std::string point_mark = "ampl 3 branch";
                                IC(point_mark);
                                //left branch is higher
                                auto ampl_left_avg = (ampl_left + ampl_central) / 2;
                                trigonometric_sample.first = ampl_left_avg;
                                auto new_loss_left = check_loss(trigonometric_sample);
                                if (new_loss_left < loss_left && new_loss_left > central_loss) {
                                    std::string point_mark = "ampl 3.1 branch";
                                    IC(point_mark);
                                    loss_left = new_loss_left;
                                    ampl_left = ampl_left_avg;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    loss_central = check_loss(trigonometric_sample);
                                }
                                else {
                                    std::string point_mark = "ampl 3.2 branch";
                                    IC(point_mark);
                                    ampl_right = (ampl_central+ampl_right) / 2;
                                    trigonometric_sample.first = ampl_right;
                                    loss_right = check_loss(trigonometric_sample);
                                    ampl_central = (ampl_right + ampl_left) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    loss_central = check_loss(trigonometric_sample);
                                }
                            }
                            else if (left_diff <= right_diff) {
                                std::string point_mark = "ampl 4 branch";
                                IC(point_mark);
                                auto ampl_right_avg = (ampl_right + ampl_central) / 2;
                                trigonometric_sample.first = ampl_right_avg;
                                auto new_loss_right = check_loss(trigonometric_sample);
                                if (new_loss_right < loss_right && new_loss_right > loss_central) {
                                    std::string point_mark = "ampl 4.1 branch";
                                    IC(point_mark);
                                    loss_right = new_loss_right;
                                    ampl_right = ampl_right_avg;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    loss_central = check_loss(trigonometric_sample);
                                }
                                else {
                                    std::string point_mark = "ampl 4.2 branch";
                                    IC(point_mark, i);
                                    ampl_left = (ampl_central + ampl_left)/2;
                                    trigonometric_sample.first = ampl_left;
                                    //IC(point_mark);
                                    loss_left = check_loss(trigonometric_sample);
                                    ampl_central = (ampl_right + ampl_left) / 2;
                                    trigonometric_sample.first = ampl_central;
                                    loss_central = check_loss(trigonometric_sample);
                                    //IC(ampl_left, ampl_central, ampl_right, loss_left, loss_central, loss_right);
                                    //IC(point_mark);
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
                                            //IC("aaaa fuck x4");
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
                            if (ampl_central - ampl_left < 0.000001 || ampl_right - ampl_central < 0.000001){
                                break;
                            }
                        }
                        if (loss_central <= loss_right && loss_central <= loss_left){
                            trigonometric_sample.first = ampl_central;
                        }
                        else if (loss_right<loss_left){
                            trigonometric_sample.first = ampl_right;
                        }
                        else{
                            trigonometric_sample.first = ampl_left;
                        }
                        
                        
                        check_loss(trigonometric_sample);
                        //show(PlottingKind::Simple);
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

