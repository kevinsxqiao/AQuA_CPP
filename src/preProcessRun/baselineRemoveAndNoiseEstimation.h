//
// Created by Kevin Qiao on 3/28/2023.
//

#ifndef AQUA_CPP_BASELINEREMOVEANDNOISEESTIMATION_H
#define AQUA_CPP_BASELINEREMOVEANDNOISEESTIMATION_H

#include "../data/data.h"


namespace AQuA{

    vector<cv::Mat> fit_F0_var(const vector<cv::Mat>& F0ProOrg, const vector<cv::Mat>& varMapOrg, int dist=0);
    vector<vector<cv::Mat>> movmean(const vector<vector<cv::Mat>>& dataIn);
    vector<cv::Mat> truncated_kept_var(const vector<cv::Mat>& quantiles);
    vector<vector<cv::Mat>> baselineLinearEstimate(vector<vector<cv::Mat>>& data);
    void correctBoundaryStd();
    float obtainBias();
    void noiseEstimationFunction(const vector<vector<cv::Mat>>& dataOrg, const vector<vector<cv::Mat>>& dataSmo,
                                 const vector<cv::Mat>& F0Pro, bool*** evtSpatialMask, vector<cv::Mat>& stdMapOrg, vector<cv::Mat>& stdMapSmo,
                                 vector<cv::Mat>& tempVarOrg,   vector<cv::Mat>& correctPars);
    vector<vector<cv::Mat>> baselineRemoveAndNoiseEstimation(vector<vector<cv::Mat>>& dataOrg, bool*** evtSpatialMask);


}// namespace

#endif //AQUA_CPP_BASELINEREMOVEANDNOISEESTIMATION_H
