module;

//#include <monster.hpp>
import modes_extractors;
import npdsp_concepts;
import <cmath>;
import <vector>;
import <numbers>;

export module encodings_generators;

using namespace NP_DSP::ONE_D;

export struct SimpleTokenizer{
    //size_t encodding_size;
    double max_value;
    double n = 10000;
    int step = 1;
    int modes_count;
    int dims_count = 8;
    std::vector<int> encode_sizes = {2, 2, 2, 10, 10, 10, 10, 10};
    std::vector<double> max_values = {500, 100, 22, 10000000000, 1.0, 1.0, 10000000000, 1570};
    bool code_metadata = true;

    template<typename T, SignalBase EncoddingT> 
    void EncodeScalar(T value, EncoddingT & encodding, size_t dim_number){
        for (int i = 0; i < encodding_size; i++){
            if (i % 2 == 0){
                encodding[i] = std::sin(value / std::pow(max_values[i]/2.0 / std::numbers::pi, i / encode_sizes[dim_number]));
            }
            else{
                encodding[i] = std::cos(value / std::pow(max_values[i]/2.0 / std::numbers::pi, (i - 1) / encode_sizes[dim_number]));
            }
        }
    }

    template<typename U, typename Conf>
    U getScalarValue(const Conf & config, size_t channel_number, size_t mode_number, size_t idx, size_t dim_couters){
        if (dim_couter == 0){
            return idx;
        } 
        else if  (dim_couter == 1){
            return mode_number;
        }
        else if  (dim_couter == 2){
            return channel_number;
        }
        else if  (dim_couter == 3){
            return (*config.modes[mode_number])[idx];
        }
        else if  (dim_couter == 4){
            return (*config.inst_freq[mode_number])[idx];
        }
        else if  (dim_couter == 5){
            return (*config.inst_ampl[mode_number])[idx];
        }
        else if  (dim_couter == 6){
            return (*config.phases[mode_number])[idx];
        }
        else if (dim_couter == -1){
            return 7;
        }
        else if (dim_couter == -2){
            return 3;
        }
    }

    template<typename U, typename Conf>
    U getScalarValueDouble(const Conf & config, size_t mode_number, size_t idx, size_t dim_couter){
        if  (dim_couter == 0){
            return idx;
        } 
        else if  (dim_couter == 1){
            return mode_number;
        } 
        else if  (dim_couter == 2){
            return channel_number;
        }
        else if  (dim_couter == 3){
            return (*config.modes[mode_number])[idx];
        }
        else if  (dim_couter == 4){
            return (*config.inst_freq[mode_number])[idx];
        }
        else if  (dim_couter == 5){
            return (*config.inst_freq[mode_number])[idx].first;
        }
        else if  (dim_couter == 6){
            return (*config.inst_ampl[mode_number])[idx].second;
        }
        else if  (dim_couter == 7){
            return (*config.phases[mode_number])[idx];
        }
        else if (dim_couter == -1){
            return 8;
        }
        else if (dim_couter == -2){
            return 3;
        }
    }

    template<typename U, typename Conf>
    U getScalarValueDoubleAsMono(const Conf & config, size_t mode_number, size_t idx, dim_couter){
        if  (dim_couter == 0){
            return idx;
        } 
        else if  (dim_couter == 1){
            return mode_number;
        } 
        else if (dim_couter == 2){
            return channel_number;
        }
        else if (dim_couter == 3){
            return (*config.modes[mode_number])[idx];
        }
        else if (dim_couter == 4){
            return 0.5 * ((*config.inst_freq[mode_number])[idx].first + (*config.inst_freq[mode_number])[idx].second);
        }
        else if (dim_couter == 5){
            return (*config.inst_ampl[mode_number])[idx];
        }
        else if (dim_couter == 6){
            return (*config.phases[mode_number])[idx];
        }
        else if (dim_couter == -1){
            return 7;
        }
        else if (dim_couter == -2){
            return 3;
        }
    }


    //OutT is SignalBase of SignalBases of Signals
    template<typename U, typename Conf, SignalBase OutT, typename Getter>
    void getModeSampleScalarWithMedata(Conf & config, OutT & out, Getter & getter, size_t idx, size_t mode_number, size_t dims){
        for (int i = 0; i < getter(-1); i++){
            if (i == dims){
                break;
            }
            out[i] = getter(config, mode_number, idx, i);
        }
    }

    template<typename U, typename Conf, SignalBase OutT, typename Getter>
    void getModeSampleScalarWithNoMedata(Conf & config, OutT & out, Getter & getter, size_t idx, size_t mode_number, size_t dims){
        for (int i = getter(-2); i < getter(-1); i++){
            if (i == dims){
                break;
            }
            out[i] = getter(config, mode_number, idx, i);
        }
    }

    template<typename U, typename Conf, SignalBase OutT, typename Getter, typename Getter2>
    void getSampleScalar(Conf & config, OutT & out, Getter & getter_ext, Getter2 & getter_int, size_t idx, size_t dims_low, size_t modes_count){
        auto size_fn_first = [&](){
            return dims_low + getter(-2);
        }
        auto ref_expr_first = [&](size_t idx){
            return out[idx];
        }
        ExpressionWrapper<U, size_t, decltype(ref_expr_first), decltype(ref_expr_first), decltype(size_expr),
            true> first_enc_expr(ref_expr_first, ref_expr_first, size_expr);
        
        auto size_fn = [&](){
            return dims_low;
        }

        getter(config, first_enc_expr, getter_int, idx, 0, dims_low + 2);
        for (int i = 1; i < modes_count; i++){
            auto ref_expr = [&](size_t idx){
                return out[idx + dims_low * i + getter(-2)];
            }
            ExpressionWrapper<U, size_t, decltype(ref_expr), decltype(ref_expr), decltype(size_expr),
            true> enc_expr(ref_expr, ref_expr, size_expr);
            getter(config, enc_expr, getter_int, idx, i, dims_low);
        }
    }

    template<SignalBase ScalarSampleT, SignalBase OutT>
    void scalarSampleEncode(ScalarSampleT & sample, OutT & out){
        auto size_expr = [&](){
            return encodding_size;
        };
        if (code_metadata){
            for(int i = 0; i < sample.size()){
                auto refExpr = [&](size_t idx){
                    return out[encodding_size * i + idx];

                };

            }
        }
        else{
            //todo
        }
    }

    template<typename Conf, typename Stop, typename ComputeFunc, Signal DataT, Signal OutT>
    void computeData(Conf & config, Stop & stop, ComputeFunc & compute, DataT & data, OutT & out){
        config.load(data);
        compute(config, stop);

        
    }
};