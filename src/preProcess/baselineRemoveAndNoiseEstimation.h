//
// Created by Kevin Qiao on 3/28/2023.
//

#ifndef AQUA_CPP_BASELINEREMOVEANDNOISEESTIMATION_H
#define AQUA_CPP_BASELINEREMOVEANDNOISEESTIMATION_H

#include "../data/data.h"
#include <iostream>


namespace AQuA{

    std::vector<std::vector<cv::Mat>> movmean(const std::vector<std::vector<cv::Mat>>& dataIn);
    void baselineLinearEstimate(std::vector<std::vector<cv::Mat>>& data);
    void baselineRemoveAndNoiseEstimation(std::vector<std::vector<cv::Mat>>& data1);


}// namespace

#endif //AQUA_CPP_BASELINEREMOVEANDNOISEESTIMATION_H
