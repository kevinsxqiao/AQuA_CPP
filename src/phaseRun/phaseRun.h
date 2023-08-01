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
    void seDetection(const vector<vector<cv::Mat>>& dF, const vector<vector<cv::Mat>>& dataOrg,
                     const vector<vector<int>>& arLst);
    vector<vector<vector<cv::Mat>>> normalizeAndResize(const vector<vector<cv::Mat>>& dataOrg);
    Score_struct getSeedScore_DS4(vector<Point_struct> pix, vector<vector<cv::Mat>> datVec, int H0, int W0, int t_scl);
    void seedDetect2_DS_accelerate(vector<vector<cv::Mat>> dF, const vector<vector<cv::Mat>>& dataOrg,
                                   const vector<vector<int>>& arLst);

}

#endif //AQUA_PHASERUN_H
