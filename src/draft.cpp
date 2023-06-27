#include <iostream>
#include <vector>
#include <algorithm>
#include "data/data.h"
//#include <boost/math/distributions/normal.hpp>

bool*** createEvt(){
    bool*** evtSpatialMask;
    evtSpatialMask = new bool** [3];
    for (int i = 0; i < 3; ++i) {
        evtSpatialMask[i] = new bool* [3];
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            evtSpatialMask[i][j] = new bool [3];
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                evtSpatialMask[i][j][k] = true;
            }
        }
    }
    return evtSpatialMask;
}

void test(bool*** a){
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                if (a[i][j][k]) {
                    std::cout << 1 << " ";
                }
            }
        }
    }
}

int main() {
    bool*** a = createEvt();
    std::cout<<1;


    return 0;
}
