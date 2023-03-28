//
// Created by Kevin Qiao on 2/17/2023.
//

#ifndef AQUA_CPP_REGCROSSCORRELATION_H
#define AQUA_CPP_REGCROSSCORRELATION_H

#include "../data/data.h"

#define N0 AQuA::rawDataSize::size1
#define N1 AQuA::rawDataSize::size2
#define N2 AQuA::rawDataSize::size3
#define N0_ext (2*N0-1)
#define N1_ext (2*N1-1)
#define N2_ext (2*N2-1)
#define FRAME AQuA::rawDataSize::frame

namespace AQuA{

    float medianFunc(float array[], int size);
    void rotate_2d(float ref[N0][N1][N2],
                   float (&mat)[N0][N1][N2]);
    void flip_3d(float ref[N0][N1][N2],
                 float (&mat)[N0][N1][N2]);
    void dft(float a_add[N0_ext][N1_ext][N2_ext],
             float b_add[N0_ext][N1_ext][N2_ext],
             float(&c)[N0_ext][N1_ext][N2_ext]);
    void dft(double a_add[N0_ext][N1_ext][N2_ext],
             double b_add[N0_ext][N1_ext][N2_ext],
             double(&c)[N0_ext][N1_ext][N2_ext]);
    void calCC(float a[N0][N1][N2],     // input: a = moving[]; b = ref[]
               float b[N0][N1][N2],     // output: c[]
               float (&c)[N0_ext][N1_ext][N2_ext]);
    float **** regCrossCorrelation(float (&data1)[N0][N1][N2][FRAME],
                                   float (&data2)[N0][N1][N2][FRAME]);




}// namespace

#endif //AQUA_CPP_REGCROSSCORRELATION_H
