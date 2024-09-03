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
}