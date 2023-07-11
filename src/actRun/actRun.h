//
// Created by Kevin Qiao on 2/10/2023.
//

#ifndef AQUA_ACTRUN_H
#define AQUA_ACTRUN_H

#include "../data/data.h"

namespace AQuA{
    void actRun();
    std::vector<std::vector<Point_struct>> bw2Reg(std::vector<std::vector<cv::Mat>> BW);
    std::vector<std::vector<Point_struct>> bwconncomp4D(std::vector<std::vector<cv::Mat>> BW);

}

#endif //AQUA_ACTRUN_H
