//
// Created by Kevin Qiao on 2/13/2023.
//

#include <iostream>
#include "data/data.h"
#include <algorithm>
#include "draft.h"
#include <fftw3.h>

#define N0 2
#define N1 3
#define N2 2
#define FRAME AQuA::rawDataSize::frame

#define REAL 0
#define IMAG 1
#define N 2

//void rotate2dMatrix(float ref[N0][N1][N2],
//                    float (&mat)[N0][N1][N2]){
//    for (int i = 0; i < N0; ++i) {
//        for (int j = 0; j < N1; ++j) {
//            mat[N0 - j][i][0] = ref[i][j][0];
//        }// for(j)
//    }// for(i)
//}// rotate2dMatrix()
//
//void flip_3d(float ref[N0][N1][N2],
//             float (&mat)[N0][N1][N2]){
//    for (int i = 0; i < N0; ++i) {
//        for (int j = 0; j < N1; ++j) {
//            for (int k = 0; k < N2; ++k) {
//                mat[N0 - 1 - i][N1 - 1 -j ][N2 - 1 - k] = ref[i][j][k];
//            }// for(k)
//        }// for(j)
//    }// for(i)
//}// flip_3d()


int main() {
    double a=2.4328753235,b=1.52435293876896,c;
    c = a/b;
    std::cout<<c<<std::endl;
//    int id;
//    float matrix[N0][N1][N2]={1,2,3,4,5,6,7,8,9,10,11,12};
//    float mat[N0][N1][N2] = {0};
////    cout<<*max_element(matrix,matrix+8)<<endl;
//    for (int i = 0; i < N0; ++i) {
//        for (int j = 0; j < N1; ++j) {
//            for (int k = 0; k < N2; ++k) {
//                std::cout<< matrix[i][j][k]<< " ";
//            }
//
//        }
//
//    }
//    std::cout<< std::endl;
//    flip_3d(matrix,mat);
//    for (int i = 0; i < N0; ++i) {
//        for (int j = 0; j < N1; ++j) {
//            for (int k = 0; k < N2; ++k) {
//                std::cout<< mat[i][j][k]<< " ";
//            }
//
//        }
//
//    }
//    std::cout<< std::endl;
//    id = max_element(matrix,matrix+8)-matrix;
//    int x;
//    x = std::min(h,w);
//    std::cout<<x;

//    cout<<id<<endl;
//    h = id/4;
//    w = id%4;
//    cout<<h<<" "<<w<<endl;

//    int a[2][3] = {1,2,3,4,5,6};
//
//    int row=0, col=0;
//    row = 4/3;
//    col = 4%3;
//    cout<<row<<endl;
//    cout<<col;
    return 0;
}



//int main() {
//
//
////    AQuA::regCrossCorrelation(AQuA::datOrg1,AQuA::datOrg2);
//
//
//
//    return 0;
//}
