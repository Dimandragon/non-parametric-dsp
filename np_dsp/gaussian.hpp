#pragma once
#include <fast_gaussian_blur_template.h>
#include <new>
#include <cstddef>

template<typename T>
void gaussianBlur(const T & in, T & out, double sigma, int passes = 5){
    int width = in.size(), height = in[0].size(), channels = 1;
    Border border = Border::kMirror;
    std::size_t size = width * height * channels;
    double * new_image = new double[size];
    double * old_image = new double[size];

    for (int i = 0; i < width; i++){
        for (int j = 0; j < height; j++){
            old_image[i * height + j] = in[i][j];
        }
    }

    fast_gaussian_blur(old_image, new_image, width, height, channels, sigma, passes, border);
    for (int i = 0; i < width; i++){
        for (int j = 0; j < height; j++){
            out[i][j] = new_image[i * height + j];
        }
    }
    delete [] new_image;
    delete [] old_image;
}