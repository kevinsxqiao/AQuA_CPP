//
// Created by Kevin Qiao on 3/28/2023.
//

#include "preProcessRun.h"

namespace AQuA{

    void preProcessRun(DATA_TYPE**** data1, DATA_TYPE**** data2){
        if(isDefault() || opts.alreadyPreprocess==0 ){ // judge whether this step is already done, since this is time-consuming
            data1 = regCrossCorrelation(data1,data2);
            /*
             *     scl.hrg = [1,size(datOrg1,1)];
                    scl.wrg = [1,size(datOrg1,2)];
                    scl.lrg = [1,size(datOrg1,3)];
             */

        }

        bool ** evtSpatialMask;
        evtSpatialMask = new bool * [H];
        for (int i = 0; i < H; ++i) {
            evtSpatialMask[i] = new bool [W];
        }
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                evtSpatialMask[i][j] = true;
            }
        }

//        % F0 bias calculation
        opts.movAvgWin = std::min(opts.movAvgWin, 100);
        opts.cut = std::min(opts.cut, 10000);
//        % load setting
//        % obtain data

//        % smooth + noise estimation + remove background
        data1 = baselineRemoveAndNoiseEstimation(data1);

    }// preProcessRun()

}// namespace

