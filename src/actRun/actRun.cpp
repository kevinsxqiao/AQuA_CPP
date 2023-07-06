//
// Created by Kevin Qiao on 7/6/2023.
//
#include "actRun.h"

namespace AQuA{

    void bw2Reg(std::vector<std::vector<cv::Mat>> BW){
        if (opts.spaMergeDist>0){
            if (BW[0].size() == 1){

            }
        }

    }


    void acDetect(std::vector<std::vector<cv::Mat>> dF1, bool*** evtSpatialMask) {
//        if (ch == 1){
//
//        }
        std::vector<float> thrs;
        if ((opts.thrARScl > opts.maxdF1) || (opts.maxSize >= H * W * L) && (opts.circularityThr == 0)) {
            //no advanced filter setting, single threshold
            thrs.push_back(opts.thrARScl);
        } else {
            //have advanced filter setting, multiple threshold
            float step = (opts.maxdF1 - opts.thrARScl) / 10;
            for (float thr = opts.thrARScl; thr <= opts.maxdF1; thr += step) {
                thrs.push_back(thr);
            }
        }//else
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    if (!evtSpatialMask[i][j][k]){
                        for (int t = 0; t < T; ++t) {
                            dF1[t][k].at<float>(i,j) = -1;
                        }//for(t)
                    }//if
                }//for(k)
            }
        }//for(i)

        //valid region
        std::vector<std::vector<cv::Mat>> activeMap(T);
        std::vector<std::vector<cv::Mat>> selectMap(T);
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                activeMap[t][k] = cv::Mat::zeros(H,W,CV_32S);
                selectMap[t][k] = cv::Mat::zeros(H,W,CV_32S);
            }
        }
        float nReg = 0;
        for (int k_thr = 0; k_thr < thrs.size(); ++k_thr) {
            float thr = thrs[k_thr];
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    for (int i = 0; i < H; ++i) {
                        for (int j = 0; j < W; ++j) {
                            if ((dF1[t][k].at<float>(i,j) > thr) && (activeMap[t][k].at<int>(i,j) == 0)){
                                selectMap[t][k].at<int>(i,j) = 1;
                            } else{
                                selectMap[t][k].at<int>(i,j) = 0;
                            }
                        }
                    }
                }
            }
        }

    }



    /*
     * active region detection and update overlay map
     */
    void actRun(){
        std::cout<< "--------start detecting--------"<<std::endl;
        bool*** evtSpatialMask = createEvtSpatialMask();

        //foreground and seed detection
        acDetect(opts.dF1, evtSpatialMask);

    }

}
