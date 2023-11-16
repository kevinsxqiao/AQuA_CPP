//
// Created by Kevin Qiao on 2/10/2023.
//

#ifndef AQUA_PHASERUN_H
#define AQUA_PHASERUN_H

#include "../data/data.h"
#include "../actRun/actRun.h"
#include "../preProcessRun/baselineRemoveAndNoiseEstimation.h"

namespace AQuA{
    void phaseRun();
    cv::Mat myResize(const cv::Mat& src, float scaleRatio_x, float scaleRatio_y);
    void seDetection(vector<vector<cv::Mat>> dF, const vector<vector<cv::Mat>>& dataOrg, vector<vector<int>>& arLst);
    vector<vector<vector<cv::Mat>>> normalizeAndResize(const vector<vector<cv::Mat>>& dataOrg);
    Score_struct getSeedScore_DS4(const vector<ushort>& pix, const vector<vector<cv::Mat>>& datVec, int H, int W, int L, int T, float t_scl);
    void seedDetect2_DS_accelerate(vector<vector<cv::Mat>> dF, const vector<vector<cv::Mat>>& dataOrg,
                                   vector<vector<int>>& arLst, vector<vector<cv::Mat>>& Map);
    void ordStatSmallSampleWith0s(const vector<float>& fg, const vector<float>& bg, const vector<float>& nanVec, double& mu, double& sigma);

}

#endif //AQUA_PHASERUN_H
