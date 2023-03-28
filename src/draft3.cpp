#include <iostream>
#include <cmath>
#include <fftw3.h>
#include<time.h>
#include <algorithm>

double *fun(double* pt) {
    for (int i = 0; i < 3; ++i) {
        pt[i]++;
    }
    return pt;
}


int main() {
    double ****pt;
    int x = 2, y = 3, z = 2;
    int index = 0;

    pt = new double ***[x];
    for (int i = 0; i < x; ++i) {
        pt[i] = new double**[y];

    }
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            pt[i][j] = new double*[z];
        }
    }
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            for (int k = 0; k < z; ++k) {
                pt[i][j][k] = new double[5];
            }
        }
    }
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            for (int k = 0; k < z; ++k) {
                for (int t = 0; t < 5; ++t) {
                    pt[i][j][k][t] = index++;
                }
            }
        }
    }
//    for (int i = 0; i < 3; ++i) {
//        pt[i] = i;
//    }
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            for (int k = 0; k < z; ++k) {
                for (int t = 0; t < 5; ++t) {
                    std::cout<<i<<","<<j<<","<<k<<":";
                    std::cout << pt[i][j][k][t] << " ";
                }
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    return 0;
}
//    fun(pt);
//    for (int i = 0; i < 3; ++i) {
//        std::cout << pt[i] << " ";
//    }
//        return 0;
//    }
//    int a[2][2][2] = {1, 2, 3, 4, 5, 6, 7, 8};
//    int b[8];
//    int index = 0;
//    int id = 0, r1;
//    int ix, iy, iz;
//    for (int i = 0; i < 2; ++i) {
//        for (int j = 0; j < 2; ++j) {
//            for (int k = 0; k < 2; ++k) {
//                b[index++] = a[i][j][k];
//            }
//        }
//    }
//    int max = std::max_element(b, b + 8) - b;
//    int min = std::min_element(b, b + 8) - b;
//    std::cout << max << std::endl;
//    std::cout << min << std::endl;
//    for (int i = 0; i < 8; ++i) {
//        std::cout << b[i] << " ";
//    }
//    std::cout << std::endl;
//    ix = id / 4;
//    r1 = id % 4;
//    iy = r1 / 2;
//    iz = r1 % 2;
//    std::cout << a[ix][iy][iz];
//    return 0;
//}

