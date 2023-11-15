//
// Created by Kevin Qiao on 2/17/2023.
//

#ifndef AQUA_CPP_REGCROSSCORRELATION_H
#define AQUA_CPP_REGCROSSCORRELATION_H

#include "../data/data.h"
#include <algorithm>
#include <fftw3.h>
#include <iostream>
#include <time.h>

namespace AQuA{

//    float medianFunc(float* array, int size);
    float*** dft(float*** a_add, float*** b_add, int H_ext, int W_ext, int L_ext);
    float*** flip3dMatrix(float*** ref, float*** b_flip, int H, int W, int L);
    float*** calCC(float*** a, float*** b, float*** a_add, float*** b_add, float*** b_flip, int H, int W, int L);
    vector<vector<cv::Mat>> regCrossCorrelation(vector<vector<cv::Mat>>& data1);
//    DATA_TYPE medianFunc(DATA_TYPE *array, int size);
//    float*** dft(float*** a_add, float*** b_add);
//    double*** dft(double*** a_add, double*** b_add);
//    DATA_TYPE*** rotate2d(DATA_TYPE*** ref);
//    DATA_TYPE*** flip3d(DATA_TYPE*** ref, DATA_TYPE*** b_flip);
//    DATA_TYPE*** calCC(DATA_TYPE*** a, DATA_TYPE*** b, DATA_TYPE*** a_add, DATA_TYPE*** b_add, DATA_TYPE*** b_flip);
//    DATA_TYPE**** regCrossCorrelation(DATA_TYPE**** data1, DATA_TYPE**** data2);

}// namespace

#endif //AQUA_CPP_REGCROSSCORRELATION_H
