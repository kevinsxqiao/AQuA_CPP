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

    DATA_TYPE medianFunc(DATA_TYPE array[], int size);
    void rotate2d(DATA_TYPE ref[N0][N1][N2],
                   DATA_TYPE (&mat)[N0][N1][N2]);
    void flip3d(DATA_TYPE ref[N0][N1][N2],
                 DATA_TYPE (&mat)[N0][N1][N2]);
    void dft(float a_add[N0_ext][N1_ext][N2_ext],
             float b_add[N0_ext][N1_ext][N2_ext],
             float(&c)[N0_ext][N1_ext][N2_ext]);
    void dft(double a_add[N0_ext][N1_ext][N2_ext],
             double b_add[N0_ext][N1_ext][N2_ext],
             double(&c)[N0_ext][N1_ext][N2_ext]);
    void calCC(DATA_TYPE a[N0][N1][N2],     // input: a = moving[]; b = ref[]
               DATA_TYPE b[N0][N1][N2],     // output: c[]
               DATA_TYPE (&c)[N0_ext][N1_ext][N2_ext]);
    DATA_TYPE **** regCrossCorrelation(DATA_TYPE (&data1)[N0][N1][N2][FRAME],
                                   DATA_TYPE (&data2)[N0][N1][N2][FRAME]);


    
}// namespace

#endif //AQUA_CPP_REGCROSSCORRELATION_H
