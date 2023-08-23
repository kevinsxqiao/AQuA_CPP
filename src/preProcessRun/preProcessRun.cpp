//
// Created by Kevin Qiao on 3/28/2023.
//

#include "preProcessRun.h"

namespace AQuA{

    void preProcessRun(vector<vector<cv::Mat>>& data1){
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
        int H = data1[0][0].rows;
        int W = data1[0][0].cols;
        int L = data1[0].size();
        int T = data1.size();
        bool*** evtSpatialMask = createEvtSpatialMask(H,W,L);

        /*
         * F0 bias calculation
         */
        opts.movAvgWin = min(opts.movAvgWin, 100);
        opts.cut = min(opts.cut, 10000);

        /*
         * smooth + noise estimation + remove background
         */
        vector<vector<cv::Mat>> dF1 = baselineRemoveAndNoiseEstimation(data1, evtSpatialMask);
//        writeDataToMatFile(dF1,"C:/Users/Kevin Qiao/Desktop/AQuA_data/dF1.mat");
//        auto start = chrono::high_resolution_clock::now();
        opts.dF1 = dF1;
//        for (int t = 0; t < T; ++t) {
//            vector<cv::Mat> frame;
//            for (int k = 0; k < L; ++k) {
//                frame.push_back(dF1[t][k].clone());
//            }
//            opts.dF1.push_back(frame);
//        }

//        auto end = chrono::high_resolution_clock::now();
//        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
//        cout << "used time: " << duration/1000 << " seconds" << endl;

        release3dMatrix_bool(evtSpatialMask,H,W);

    }// preProcessRun()

}// namespace

