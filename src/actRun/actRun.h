//
// Created by Kevin Qiao on 2/10/2023.
//

#ifndef AQUA_ACTRUN_H
#define AQUA_ACTRUN_H

#include "../data/data.h"

namespace AQuA{
    void actRun();
    vector<vector<int>> bw2Reg(const vector<vector<cv::Mat>>& BW);
    vector<vector<int>> bwconncomp4D(const vector<vector<cv::Mat>>& BW);
    vector<vector<Point_struct>> acDetect(vector<vector<cv::Mat>>& dF1, bool*** evtSpatialMask);

}

#endif //AQUA_ACTRUN_H
