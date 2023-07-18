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

        bool*** evtSpatialMask = createEvtSpatialMask();

        /*
         * F0 bias calculation
         */
        opts.movAvgWin = std::min(opts.movAvgWin, 100);
        opts.cut = std::min(opts.cut, 10000);

        /*
         * smooth + noise estimation + remove background
         */
        std::vector<std::vector<cv::Mat>> dF1 = baselineRemoveAndNoiseEstimation(data1, evtSpatialMask);
//        writeDataToMatFile(dF1,"C:/Users/Kevin Qiao/Desktop/AQuA_data/dF1.mat");
//        auto start = std::chrono::high_resolution_clock::now();
        opts.dF1 = dF1;
//        for (int t = 0; t < T; ++t) {
//            std::vector<cv::Mat> frame;
//            for (int k = 0; k < L; ++k) {
//                frame.push_back(dF1[t][k].clone());
//            }
//            opts.dF1.push_back(frame);
//        }

//        auto end = std::chrono::high_resolution_clock::now();
//        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//        std::cout << "used time: " << duration/1000 << " seconds" << std::endl;

        release3dMatrix_bool(evtSpatialMask,H,W);

    }// preProcessRun()

}// namespace

