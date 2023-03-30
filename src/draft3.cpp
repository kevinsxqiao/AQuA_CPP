#include <iostream>
#include <cmath>
#include <fftw3.h>
#include<time.h>
#include <algorithm>
#include "preProcess/regCrossCorrelation.h"



int main() {
    int a[3] = {1,2,3};
    int* b = new int[3];
    for (int i = 0,m=0; i < 3; ++i) {
        b[i]=i;
    }
    std::cout<<std::max_element(b,b+3)-b;
    return 0;
}

