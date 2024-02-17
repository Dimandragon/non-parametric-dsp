module;

export module modes_extractors;

import npdsp_concepts;

namespace NP_DSP{
    namespace ONE_D{
        namespace MODES_EXTRACTORS{
            export
            template <InstFreqComputer InstFreqComputerT>
            struct LinearExtractor{
                using InstFreqComputerT 
            }

            export 
            struct MainExtractor{
                using DataType = NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true>;
                DataType data;
                std::vector<DataTypes *> modes;
                std::vector<DataTypes *> inst_freqs;
                std::vector<DataTypes *> inst_ampls;
                std::vector<DataTypes *> phases;

                ~MainExtractor(){
                    for (int i = 0; i < modes.size(); i++){
                        delete modes[i];
                    }
                    for (int i = 0; i < inst_freqs.size(); i++){
                        delete inst_freqs[i];
                    }
                    for (int i = 0; i < phases.size(); i++){
                        delete phases[i];
                    }
                    for (int i = 0; i < inst_ampls.size(); i++){
                        delete inst_ampls[i];
                    }
                }

                template<typename DataInT>
                void load(const & DataInT data_in){
                    data.base->vec->clear();
                    for (auto i = 0; i < data_in.size(); i++){
                        data.push_back(data.size());
                    }
                    for (int i = 0; i < modes.size(); i++){
                        modes[i]->base->vec->clear();
                    }
                    for (int i = 0; i < inst_freqs.size(); i++){
                        inst_freqs[i]->base->vec->clear();
                    }
                    for (int i = 0; i < phases.size(); i++){
                        phases[i]->base->vec->clear();
                    }
                    for (int i = 0; i < inst_ampls.size(); i++){
                        inst_ampls[i]->base->vec->clear();
                    }
                }

                bool stop(size_t iter_number){
                    if ((*(phases[iter_number]))[data.size()-1] < std::numbers::pi * 2){
                        return true;
                    }
                    else{
                        return false;
                    }
                }

                void compute(){
                    size_t iter_number = 0;
                    for(;;){
                        //todo

                        if (stop(iter_number)){
                            break;
                        }
                        iter_number++;
                    }
                }
            }
        }
    }
}