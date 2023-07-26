//
// Created by Kevin Qiao on 2/10/2023.
//

#ifndef AQUA_PHASERUN_H
#define AQUA_PHASERUN_H

#include "../data/data.h"
#include "../actRun/actRun.h"

namespace AQuA{
    void phaseRun();
    cv::Mat myResize(const cv::Mat& src, float scaleRatio_x, float scaleRatio_y);

}

#endif //AQUA_PHASERUN_H
