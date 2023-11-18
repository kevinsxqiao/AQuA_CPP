#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <mat.h>
#include <vector>
#include <boost/math/distributions/normal.hpp>
#include "data/data.h"
#include <omp.h>
#include "phaseRun/phaseRun.h"
#include "actRun/actRun.h"
#include <chrono>





int main() {
    int H=190,W=190,L=50,T=200;
    vector<vector<cv::Mat>> zscoreMap_selected = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/zScore.mat", "zScore");
    vector<vector<cv::Mat>> zscoreMap = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/zscoreMap.mat", "zscoreMap");
    vector<vector<cv::Mat>> activeMap = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/activeMap.mat", "activeMap");
//    vector<vector<cv::Mat>> sel_cpp = AQuA::load4D_8U("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/sel_cpp.mat", "myVar");
//    vector<vector<cv::Mat>> dF1 = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/phaseRun.mat", "dF1");
//    AQuA::writeDataToMatFile(dF1, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test/dF_test.mat");
//    vector<vector<cv::Mat>> dF1_test = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/dF_test.mat", "myVar");
//    vector<vector<cv::Mat>> dF_cpp = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/dFRe0.mat", "myVar");
//    vector<vector<cv::Mat>> val_cpp = AQuA::load4D_8U("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/val0.mat", "myVar");
//    vector<vector<cv::Mat>> dF_mat = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/dFRe0_mat.mat", "mat11");
//    vector<vector<cv::Mat>> val_mat = AQuA::load4D_8U("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/val0_mat.mat", "mat1");
    vector<vector<cv::Mat>> Map(T,vector<cv::Mat>(L));
    float scaleRatio =8;
    for (int t = 0; t < T; ++t) {
        for (int k = 0; k < L; ++k) {
            Map[t][k] = cv::Mat::zeros(H,W,CV_16U);
        }
    }
//    vector<vector<int>> curRegions = AQuA::bw2Reg(sel_mat);

    int cnt = 0;
    vector<vector<int>> sdLst = AQuA::bw2Reg(zscoreMap_selected);
    for (int i = 0; i < sdLst.size(); ++i) {
        vector<int> pix;
        pix.reserve(sdLst[i].size());
        for (auto& elem:sdLst[i]) {//sd:st[19] different  potential cause is zscoreMap due to z_socre_min
            pix.push_back(elem);
        }
        unordered_set<float> scores;
        for (int pix_ind = 0; pix_ind < pix.size(); ++pix_ind) {
            AQuA::Point_struct pix_temp = AQuA::ind2sub(pix[pix_ind],H,W,L);
            scores.insert(zscoreMap[pix_temp.t][pix_temp.k].at<float>(pix_temp.i,pix_temp.j));
        }
        if (scores.size()>1){
            pix.clear();
            float scores_max = *max_element(scores.begin(),scores.end());
            for (auto& pix_ele:pix) {
                AQuA::Point_struct pix_ele_temp = AQuA::ind2sub(pix_ele,H,W,L);
                if (zscoreMap[pix_ele_temp.t][pix_ele_temp.k].at<float>(pix_ele_temp.i, pix_ele_temp.j) == scores_max){
                    pix.push_back(pix_ele);
                }
            }
            cnt++;
        }// if (scores.size()>1)
        for(auto& pix_ele:pix){
            AQuA::Point_struct pix_ele_temp = AQuA::ind2sub(pix_ele,H,W,L);
            Map[pix_ele_temp.t][pix_ele_temp.k].at<ushort>(pix_ele_temp.i, pix_ele_temp.j) = i+1;
        }

    }//for i = 1:numel(sdLst)
    vector<vector<int>> arLst;
    vector<vector<cv::Mat>> arLst_selected(T,vector<cv::Mat>(L));
    for (int t = 0; t < T; ++t) {
        for (int k = 0; k < L; ++k) {
            arLst_selected[t][k] = cv::Mat::zeros(H,W,CV_8U);
            cv::Mat mask = (Map[t][k] > 0) | (activeMap[t][k] >0);
            arLst_selected[t][k].setTo(1,mask);
        }
    }
    arLst = AQuA::bw2Reg(arLst_selected);
    return 0;
}

