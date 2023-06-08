//
// Created by Kevin Qiao on 3/28/2023.
//

#include "preProcessRun.h"

namespace AQuA{

    void preProcessRun(std::vector<std::vector<cv::Mat>> data1){
        if(isDefault() || !opts.alreadyPreprocess || !opts.alreadyBleachCorrect){ // judge whether this step is already done, since this is time-consuming
            /*
             * image registration
             */
            if (opts.registrateCorrect == 2){
                data1 = regCrossCorrelation(data1);
            }

            /*
             * bleach correction
             */
            if(opts.bleachCorrect == 2){

            }

            /*
             * median filter to remove salt and pepper noise
             */
            if(opts.medSmo > 0){

            }

//            opts.alreadyPreprocess = true;

        }

        /*
         * only consider the pixels in the drawn cells
         */

//        bool ** evtSpatialMask;
//        evtSpatialMask = new bool * [H];
//        for (int i = 0; i < H; ++i) {
//            evtSpatialMask[i] = new bool [W];
//        }
//        for (int i = 0; i < H; ++i) {
//            for (int j = 0; j < W; ++j) {
//                evtSpatialMask[i][j] = true;
//            }
//        }

        /*
         * F0 bias calculation
         */
        opts.movAvgWin = std::min(opts.movAvgWin, 100);
        opts.cut = std::min(opts.cut, 10000);

        /*
         * smooth + noise estimation + remove background
         */
//        data1 = baselineRemoveAndNoiseEstimation(data1);

    }// preProcessRun()

}// namespace

