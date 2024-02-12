
module; 

import approximators;
import npdsp_concepts;
import npdsp_config;
import <utility>;
import <vector>;
import <cmath>;
import <numbers>;
import <complex>;


export module phase_computers;

namespace NP_DSP{
    namespace ONE_D{
        namespace PHASE_COMPUTERS{
            export
            template <Signal DataT, Signal OutT, Signal AdditionalDataT,
                Integrator IntegratorT, Derivator DerivatorT>
            struct ExtremumsBased{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = AdditionalDataT;
                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;
                int tile_size = 128;
                double approx_order_coeff = 1.0;

                constexpr static bool is_phase_computer = true;

                ExtremumsBased(IntegratorType integrator_in, DerivatorType derivator_in){
                    integrator = integrator_in;
                    derivator = derivator_in;
                }

                void compute(const DataType & data, OutType & out, AdditionalDataT & computer_buffer)
                {
                    //compute extremums
                    //compute "support" vector of pairs size_t and double, 
                    //where first is extremum position and second is value = idx * std::numbers::pi
                    //compute f(x) = integral(|arctg'(data')|)
                    //compute a(x) -> min for each pair in support 
                    // reduce f(position) * a(position) - support)^2 
                    // + for each x reduce
                    //if f(x) * a(x) > f(x-1) * a(x-1) then (f(x) * a(x) - f(x-1) * a(x-1))^2 else 0

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


                    std::vector<std::pair<size_t, float>> support;
                    for (auto i = 0; i < extremums.size(); i++){
                        support.push_back({extremums[i], i * std::numbers::pi});
                    }

                    GENERAL::Nil nil;

                    derivator.compute(data, computer_buffer, nil);
                    for (int i = 0; i < data.size(); i++){
                        computer_buffer[i] = std::atan(computer_buffer[i]);
                    }
                    derivator.compute(computer_buffer, out, nil);
                    for (int i = 0; i < data.size(); i++){
                        out[i] = std::abs(computer_buffer[i]);
                    }
                    integrator.compute(out, computer_buffer, nil);
                    computer_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);

                    auto loss = [&](auto & approximator){
                        double error = 0.0;
                        for (int i = 0; i < support.size(); i++){
                            auto position = support[i].first;
                            auto val = support[i].second;
                            error += 
                                (val - approximator.compute(position) * computer_buffer[position])
                                * (val - approximator.compute(position) * computer_buffer[position]) 
                                * computer_buffer.size() / support.size() / computer_buffer[position];
                        }
                        for (int i = 1; i < computer_buffer.size(); i++){
                            for (int j = 1; j <= i; j++){
                                if (computer_buffer[i]*approximator.compute(i) <
                                        computer_buffer[i-j]*approximator.compute(i-j)){
                                    error += (computer_buffer[i]*approximator.compute(i) - computer_buffer[i-j]*approximator.compute(i-j)) *
                                        (computer_buffer[i]*approximator.compute(i) - computer_buffer[i-j]*approximator.compute(i-j)) / computer_buffer[i];
                                }
                            }
                        }
                        return error;
                    };

                    auto bySampleError = [&](auto & approximator, auto idx){
                        /*double error = 0.0;
                        auto i_min = 0;
                        auto i_max = extremums.size() - 1;
                        
                        while (i_max - i_min > 1){
                            auto central = i_min + i_max / 2;
                            if (extremums[central] > idx){
                                i_min = central;
                            }
                            else{
                                i_max = central;
                            }
                        }
                        auto val_left = support[i_min].second;
                        auto val_right = support[i_max].second;
                        auto id_left = support[i_min].first;
                        auto id_right = support[i_max].first;
                        error += (val_left - approximator.approximated_data[id_left].real() * computer_buffer[id_left]) 
                            * (val_left - approximator.approximated_data[id_left].real() * computer_buffer[id_left]);
                        error += (val_right - approximator.approximated_data[id_right].real() * computer_buffer[id_right]) 
                            * (val_right - approximator.approximated_data[id_right].real() * computer_buffer[id_right]);
                        error = error / 2;
                        
                        if (idx > 0){
                            if (approximator.approximated_data[idx - 1].real() * computer_buffer[idx - 1] -
                                approximator.approximated_data[idx].real() * computer_buffer[idx] > 0)
                            {
                                error += (approximator.approximated_data[idx - 1].real() * computer_buffer[idx - 1] -
                                    approximator.approximated_data[idx].real() * computer_buffer[idx]) 
                                    * (approximator.approximated_data[idx - 1].real() * computer_buffer[idx - 1] -
                                    approximator.approximated_data[idx].real() * computer_buffer[idx]);
                            }
                        }
                        if (idx < data.size() - 1){
                            if (approximator.approximated_data[idx].real() * computer_buffer[idx] -
                                approximator.approximated_data[idx + 1].real() * computer_buffer[idx + 1] > 0)
                            {
                                error += (approximator.approximated_data[idx + 1].real() * computer_buffer[idx + 1] -
                                    approximator.approximated_data[idx].real() * computer_buffer[idx]) 
                                    * (approximator.approximated_data[idx + 1].real() * computer_buffer[idx + 1] -
                                    approximator.approximated_data[idx].real() * computer_buffer[idx]);
                            }
                        }
                        return error;
                        */
                       double error = 0.0;
                        for (int i = 0; i < support.size(); i++){
                            auto position = support[i].first;
                            auto val = support[i].second;
                            error += 
                                (val - approximator.compute(position) * computer_buffer[position])
                                * (val - approximator.compute(position) * computer_buffer[position]) 
                                * computer_buffer.size() / support.size() / computer_buffer[position];
                        }
                        for (int i = 1; i < computer_buffer.size(); i++){
                            for (int j = 1; j <= i; j++){
                                if (computer_buffer[i]*approximator.compute(i) <
                                        computer_buffer[i-j]*approximator.compute(i-j)){
                                    error += (computer_buffer[i]*approximator.compute(i) - computer_buffer[i-j]*approximator.compute(i-j)) *
                                        (computer_buffer[i]*approximator.compute(i) - computer_buffer[i-j]*approximator.compute(i-j)) / computer_buffer[i];
                                }
                            }
                        }
                        return error;
                    };

                    auto stopPoint = [](auto losses_different, auto & approximator) {
                        if (losses_different > 0.00001){
                            return false;
                        }
                        else{
                            return true;
                        }
                    };

                    auto approximator = APPROX::FourierSeriesBased<DataType, decltype(loss), decltype(stopPoint), 
                        APPROX::FSApproxKind::Positive, decltype(bySampleError)>
                            (loss, out, stopPoint);
                    approximator.tile_size = tile_size;
                    approximator.bySampleLoss = &bySampleError;       
                    approximator.is_actual = false;
                    approximator.setApproxOrderRatio(approx_order_coeff);
                    approximator.train();

                    approximator.is_actual = false;
                    for (auto i = 0; i < out.size(); i++){
                        out[i] = computer_buffer[i] * approximator.compute(i);
                    }
                    //todo test it and write example
                }
            };
        }
    }
}