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



<<<<<<< HEAD
        // compute the mean value from frame_start to frame_end
        int frame_start = 0;
        int frame_end = 2;
        for (int i=0, m=0; i<N0;++i) {
            for (int j=0; j<N1;++j) {
                for(int k=0; k<N2;++k){
                    sum = 0;
                    for (int t=frame_start;t<frame_end;++t) {
                        sum += data1[i][j][k][t];
                    }// for(t)
                    ref[i][j][k] = sum / float(frame_end - frame_start);
                    temp_col[m] = ref[i][j][k]; //convert ref to one dimension
                    ++m;
                }// for(k)
            }// for(j)
        }// for(i)

        // compute the median value of temp_col --- medianFunc(<array>, <size>)
        median = medianFunc(temp_col, N0 * N1 * N2);


        // align dark
        for(int i=0;i<N0;++i){
            for(int j=0;j<N1;++j){
                for(int k=0;k<N2;++k){
                    ref[i][j][k] -= median;  // ref = ref - median(ref(:));
                }// for(k)
            }// for(j)
        }// for(i)


        // concise cross correlation --loop from frame <2> to <FRAME>
        for(int t=0; t<FRAME; ++t){
            for(int i=0, m=0; i<N0; ++i){
                for(int j=0; j<N1; ++j){
                    for(int k=0; k<N2; ++k){
                        moving[i][j][k] = data1[i][j][k][t]; //initialize moving
                        temp_col[m] = moving[i][j][k]; //convert moving to one dimension
                        ++m;
                    }// for(k)
                } // for(j)
            }// for(i)
            median = medianFunc(temp_col, N0 * N1 * N2);
            for(int i=0;i<N0;++i){
                for(int j=0;j<N1;++j){
                    for(int k=0;k<N2;++k){
                        moving[i][j][k] -= median;  // moving = moving - median(moving(:));
                        calCC(moving, ref, matrix); // matrix = calCC(moving,ref);
                    }// for(k)
                } //for(j)
            }// for(i)

            for(int i=0, m=0;i<N0_ext;++i){
                for(int j=0;j<N1_ext;++j){
                    for(int k=0;k<N2_ext;++k){
                        matrix_col[m] = matrix[i][j][k];  // convert matrix to one-dimension to find the position of the maximum element
                    }// for(k)
                } //for(j)
            }// for(i)

            int id = std::max_element(matrix_col, matrix_col+(N0_ext*N1_ext*N2_ext)) - matrix_col;  //[~,id] = max(matrix(:));
            // [hShift,wShift,lShift] = ind2sub(size(matrix),id);
            hShift = id/(N1_ext*N2_ext); // x as page, in which page
            int r1 = id%(N1_ext*N2_ext); // in the page hShift(2d matrix), linear position
            wShift = r1 / N2_ext; // in the page hShift(2d matrix), which row;
            lShift = r1 % N2_ext; // in the page hShift(2d matrix), which column;
            x_translation[t] = N0 - hShift;
            y_translation[t] = N1 - wShift;
            z_translation[t] = N2 - lShift;
        }// for(t)


        for(int t=0;t<FRAME;++t){
            if (x_translation[t]>=0){
                xs0 = 1;
                xe0 = N0 - x_translation[t];
                xs1 = 1 + x_translation[t];
                xe1 = N0;
            }
            else{
                xs0 = 1 - x_translation[t];
                xe0 = N0;
                xs1 = 1;
                xe1 = N0 + x_translation[t];
            }

            if (y_translation[t]>=0){
                ys0 = 1;
                ye0 = N1 - y_translation[t];
                ys1 = 1 + y_translation[t];
                ye1 = N1;
            }
            else{
                ys0 = 1 - y_translation[t];
                ye0 = N1;
                ys1 = 1;
                ye1 = N1 + y_translation[t];
            }

            if (z_translation[t]>=0){
                zs0 = 1;
                ze0 = N2 - z_translation[t];
                zs1 = 1 + z_translation[t];
                ze1 = N2;
            }
            else{
                zs0 = 1 - z_translation[t];
                ze0 = N2;
                zs1 = 1;
                ze1 = N2 + z_translation[t];
            }

            //data1(xs1:xe1,ys1:ye1,zs1:ze1,t) = data1(xs0:xe0,ys0:ye0,zs0:ze0,t);
            for (int i = xs0, x_target = xs1; i <= xe0; ++i, ++x_target) {
                for (int j = ys0, y_target = ys1; j <= ye0; ++j, ++y_target) {
                    for (int k = zs0, z_target = zs1; k <= zs1; ++k, ++z_target) {
                        data1[x_target][y_target][z_target][t] = data1[i][j][k][t];
                    } //for(k)
                } //for(j)
            } //for(i)
//            if(~isempty(data2))
//                data2(xs1:xe1,ys1:ye1,zs1:ze1,t) = data2(xs0:xe0,ys0:ye0,zs0:ze0,t);
//            end
        } //for(t)
        int x_max = std::max_element(x_translation, x_translation+FRAME) - x_translation;
        int x_min = std::min_element(x_translation, x_translation+FRAME) - x_translation;
        int y_max = std::max_element(y_translation, y_translation+FRAME) - y_translation;
        int y_min = std::min_element(y_translation, y_translation+FRAME) - y_translation;
        int z_max = std::max_element(z_translation, z_translation+FRAME) - z_translation;
        int z_min = std::min_element(z_translation, z_translation+FRAME) - z_translation;

        double ****dat1;
        int x = x_max - x_min + 1, y = y_max - y_min + 1, z = z_max - z_min+ 1; // x,y,z indicates the size of the new matrix after registration
        int index = 0;

        dat1 = new double ***[x];
        for (int i = 0; i < x; ++i) {
            dat1[i] = new double**[y];

        }
        for (int i = 0; i < x; ++i) {
            for (int j = 0; j < y; ++j) {
                dat1[i][j] = new double*[z];
            }
        }
        for (int i = 0; i < x; ++i) {
            for (int j = 0; j < y; ++j) {
                for (int k = 0; k < z; ++k) {
                    dat1[i][j][k] = new double[FRAME];
                }
            }
        }
//        data1 = data1(max(x_translation)+1:end+min(x_translation),max(y_translation)+1:end+min(y_translation),max(z_translation)+1:end+min(z_translation),:);
        for (int i = x_min, x_target=0; i < x_max+1; ++i, ++x_target) {
            for (int j = y_min, y_target=0; j < y_max+1; ++j, ++y_target) {
                for (int k = z_min, z_target=0; k < z_max+1; ++k, ++z_target) {
                    for (int t = 0; t < FRAME; ++t) {
                        dat1[x_target][y_target][z_target][t] = data1[i][j][k][t];
                    }
                }
            }
        }
//        if(~isempty(data2))
//            data2 = data2(max(x_translation)+1:end+min(x_translation),max(y_translation)+1:end+min(y_translation),max(z_translation)+1:end+min(z_translation),:);
//        end

    }// regCrossCorrelation()

    //1111

=======
>>>>>>> main

}// namespace

#endif //AQUA_CPP_REGCROSSCORRELATION_H
