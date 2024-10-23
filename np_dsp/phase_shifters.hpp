#pragma once

#include <utility_math.hpp>
#include <vector>
#include <npdsp_concepts.hpp>
#include <utility>

namespace NP_DSP::ONE_D::PHASE_SHIFTERS{
    struct HTBased{
        double phase_shift = 0.5 * std::numbers::pi;

        template<typename DataT, typename OutT>
        void compute(const DataT & data, OutT & out){
            std::vector<std::complex<double>> spectre;
            std::vector<std::complex<double>> signal;

            for (int i = 0; i < data.size(); i++){
                signal.push_back({data[i], 0.0});
                spectre.push_back({0.0, 0.0});
            }

            UTILITY_MATH::fftc2c<double>(signal, spectre);
            for (int i = 1; i < data.size(); i++){
                std::pair<double, double> t_sample = UTILITY_MATH::convertFSampleC2T(spectre[i]);
                double theta = t_sample.second;
                double ampl = t_sample.first;
                theta = theta + phase_shift;
                if (theta > std::numbers::pi){
                    theta = theta - std::numbers::pi;
                }
                if (theta < -std::numbers::pi){
                    theta = theta + std::numbers::pi;
                }
                spectre[i] = UTILITY_MATH::convertFSampleT2C<double>({ampl, theta});
            }
            UTILITY_MATH::ifftc2c<double>(spectre, signal);
            for (int i = 0; i < data.size(); i++){
                //out[i] = std::sqrt(signal[i].real()*signal[i].real() + signal[i].imag()*signal[i].imag());
                out[i] = signal[i].real();
            }
        }
    };

    /*struct FTBased{
        double phase_shift = 0.5 * std::numbers::pi;

        template<typename DataT, typename OutT>
        void compute(const DataT & data, OutT & out){
            std::vector<double> spectre;
            std::vector<double> signal;

            for (int i = 0; i < data.size(); i++){
                signal.push_back({data[i], 0.0});
                spectre.push_back({0.0, 0.0});
            }

            UTILITY_MATH::fftc2c<double>(signal, spectre);
            for (int i = 0; i < data.size(); i++){
                double imag = std::numbers::pi * i / (double)data.size() * 2.0;
                std::complex<double> muller = (0.0, std::numbers::pi * i / data.size());
                muller = std::pow(muller, power);
                spectre[i] = spectre[i] * muller;
            }
            UTILITY_MATH::ifftc2c<double>(spectre, signal);
        }
    };*/

    struct NaiveExtremumsPhaseShifter{
        double phase_shift = 0.5 * std::numbers::pi;

        template<typename DataT, typename OutT, typename PhaseT>
        void compute(const DataT & data, OutT & out, PhaseT & phase){

        }
    };
}
