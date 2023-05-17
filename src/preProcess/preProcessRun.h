//
// Created by Kevin Qiao on 2/10/2023.
//

#ifndef AQUA_PREPROCESSRUN_CPP
#define AQUA_PREPROCESSRUN_CPP

#include "../data/data.h"
#include "regCrossCorrelation.h"
#include "baselineRemoveAndNoiseEstimation.h"
#include <cstdlib>
#include <iostream>

namespace AQuA{

    void preProcessRun(float**** data1, float**** data2);

}// namespace

#endif //AQUA_PREPROCESSRUN_CPP
