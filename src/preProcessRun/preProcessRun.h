//
// Created by Kevin Qiao on 2/10/2023.
//

#ifndef AQUA_PREPROCESSRUN_CPP
#define AQUA_PREPROCESSRUN_CPP

#include "../data/data.h"
#include "regCrossCorrelation.h"
#include "baselineRemoveAndNoiseEstimation.h"
#include <cstdlib>


namespace AQuA{

    std::vector<std::vector<cv::Mat>> preProcessRun(std::vector<std::vector<cv::Mat>> data1);

}// namespace

#endif //AQUA_PREPROCESSRUN_CPP
