//
// Created by Kevin Qiao on 3/28/2023.
//

#ifndef AQUA_CPP_BASELINEREMOVEANDNOISEESTIMATION_H
#define AQUA_CPP_BASELINEREMOVEANDNOISEESTIMATION_H

#include "../data/data.h"


namespace AQuA{

    bool*** createEvtSpatialMask();
    std::vector<cv::Mat> fit_F0_var(const std::vector<cv::Mat>& F0ProOrg, const std::vector<cv::Mat>& varMapOrg);
    std::vector<std::vector<cv::Mat>> movmean(const std::vector<std::vector<cv::Mat>>& dataIn);
    std::vector<cv::Mat> truncated_kept_var(const std::vector<cv::Mat>& quantiles);
    std::vector<std::vector<cv::Mat>> baselineLinearEstimate(std::vector<std::vector<cv::Mat>>& data);
    void correctBoundaryStd();
    float obtainBias();
    void noiseEstimationFunction(const std::vector<std::vector<cv::Mat>>& dataOrg, const std::vector<std::vector<cv::Mat>>& dataSmo,
                                 const std::vector<cv::Mat>& F0Pro, bool*** evtSpatialMask, std::vector<cv::Mat>& stdMapOrg, std::vector<cv::Mat>& stdMapSmo,
                                 std::vector<cv::Mat>& tempVarOrg,   std::vector<cv::Mat>& correctPars);
    std::vector<std::vector<cv::Mat>> baselineRemoveAndNoiseEstimation(std::vector<std::vector<cv::Mat>>& dataOrg, bool*** evtSpatialMask);


}// namespace

#endif //AQUA_CPP_BASELINEREMOVEANDNOISEESTIMATION_H
