#pragma once
#include <npdsp_concepts.hpp>
#include <math.h>
#include <icecream.hpp>

template<typename Data1T, typename Data2T>
double computeOrthogonality(const Data1T & data1, const Data2T & data2){
    double avg1 = 0.0;
    double avg2 = 0.0;

    for (int i = 0; i < data1.size(); i++){
        avg1 += data1[i] / data1.size();
        avg2 += data2[i] / data2.size();
    }

    auto data1_normed = [&](auto idx){
        return data1[idx] - avg1;
    };
    auto data2_normed = [&](auto idx){
        return data2[idx] - avg2;
    };

    double data1_l2 = 0.0;
    double data2_l2 = 0.0;
    //normed as /size()
    for (auto i = 0; i < data1.size(); i++){
        data1_l2 += data1_normed(i)*data1_normed(i)/(double)(data1.size() * data1.size());
        data2_l2 += data2_normed(i)*data2_normed(i)/(double)(data2.size()*data2.size());
    }

    data1_l2 = std::sqrt(data1_l2);
    data2_l2 = std::sqrt(data2_l2);

    double normed_dot = 0.0;
    for (auto i = 0; i < data1.size(); i++){
        normed_dot += data1_normed(i) * data2_normed(i)/data1_l2 /data2_l2/data1.size()/data2.size();
    }
    //IC(data1_l2, data2_l2, data1.size(), data2.size());

    return normed_dot;
}/*for(int i = 0; i < extractor.modes.size(); i++) {
        std::stringstream path;
        path << "/home/dmitry/projects/non-parametric-dsp/examples/modes_extractors/mode" << i << ".svg";

        if (save){
            extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple, path.str());
        }
        else{
            extractor.modes[i]->show(NP_DSP::ONE_D::PlottingKind::Simple);
        }
    }*/

    /*
    std::vector<std::vector<double>> orthogonality_distr;


    std::cout << computeOrthogonality<decltype(data), decltype(data)>
            (*(extractor.modes[0]), *(extractor.modes[1])) << " " 
            << computeOrthogonality<decltype(data), decltype(data)>
            (*(extractor.modes[1]), *(extractor.modes[0])) << std::endl;

    for (int i = 0; i < extractor.modes.size(); i++){
        orthogonality_distr.push_back({});
        for (int j = 0; j < extractor.modes.size(); j++){
            //std::cout << i << " " << j << std::endl;

            std::cout << computeOrthogonality<decltype(data), decltype(data)>
            (*(extractor.modes[i]), *(extractor.modes[j])) << " ";
            
            orthogonality_distr[i].push_back(std::round(100. * computeOrthogonality<decltype(data), decltype(data)>
            (*(extractor.modes[i]), *(extractor.modes[j]))) / 100.);
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    matplot::heatmap(orthogonality_distr);
    matplot::show();
    if (save){
        matplot::save("/home/dmitry/projects/non-parametric-dsp/examples/modes_extractors/orthogonality_test.svg");
    }

    */