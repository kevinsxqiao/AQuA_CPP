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

    DATA_TYPE medianFunc(DATA_TYPE *array, int size);
    float*** dft(float*** a_add, float*** b_add);
    double*** dft(double*** a_add, double*** b_add);
    DATA_TYPE*** rotate2d(DATA_TYPE*** ref);
    DATA_TYPE*** flip3d(DATA_TYPE*** ref, DATA_TYPE*** b_flip);
    DATA_TYPE*** calCC(DATA_TYPE*** a, DATA_TYPE*** b, DATA_TYPE*** a_add, DATA_TYPE*** b_add, DATA_TYPE*** b_flip);
    DATA_TYPE**** regCrossCorrelation(DATA_TYPE**** data1, DATA_TYPE**** data2);

}// namespace

#endif //AQUA_CPP_REGCROSSCORRELATION_H
