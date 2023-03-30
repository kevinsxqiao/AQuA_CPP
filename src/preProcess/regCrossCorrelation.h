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
    void rotate2d(DATA_TYPE ref[H][W][L],
                   DATA_TYPE (&mat)[H][W][L]);
    void flip3d(DATA_TYPE ref[H][W][L],
                 DATA_TYPE (&mat)[H][W][L]);
    void dft(float a_add[H_ext][W_ext][L_ext],
             float b_add[H_ext][W_ext][L_ext],
             float(&c)[H_ext][W_ext][L_ext]);
    void dft(double a_add[H_ext][W_ext][L_ext],
             double b_add[H_ext][W_ext][L_ext],
             double(&c)[H_ext][W_ext][L_ext]);
    void calCC(DATA_TYPE a[H][W][L],     // input: a = moving[]; b = ref[]
               DATA_TYPE b[H][W][L],     // output: c[]
               DATA_TYPE (&c)[H_ext][W_ext][L_ext]);
    DATA_TYPE **** regCrossCorrelation(DATA_TYPE (****data1),
                                   DATA_TYPE (****data2));


    
}// namespace

#endif //AQUA_CPP_REGCROSSCORRELATION_H
