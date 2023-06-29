#include <iostream>
#include <vector>
#include <algorithm>
#include "data/data.h"
#include "preProcess/preProcessRun.h"



int main() {
    AQuA::Init();
    float bias = AQuA::obtainBias();
    std::cout<<bias;
    return 0;
    }




