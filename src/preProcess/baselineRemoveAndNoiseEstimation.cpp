//
// Created by Kevin Qiao on 3/28/2023.
//

#include "baselineRemoveAndNoiseEstimation.h"


namespace AQuA{

    std::vector<std::vector<cv::Mat>> movmean(const std::vector<std::vector<cv::Mat>>& dataIn) {
        std::vector<std::vector<cv::Mat>> dataOut(T, std::vector<cv::Mat>(L));
        int halfWin = opts.movAvgWin/2;
        for (int k = 0; k < L; ++k) {
            for (int t = 0; t < T; ++t) {
                int start = std::max(0, t - halfWin);
                int end = std::min(T - 1, t + halfWin);
                int count = 0;
                cv::Mat sum = cv::Mat::zeros(dataIn[0][0].rows, dataIn[0][0].cols, CV_32F);
                for (int i = start; i <= end; ++i) {
                    cv::add(sum, dataIn[i][k], sum, cv::noArray(), -1);
                    ++count;
                }
                dataOut[t][k] = sum / count;
            }
        }
        return dataOut;
    }


    void baselineLinearEstimate(std::vector<std::vector<cv::Mat>>& data){
        std::vector<std::vector<cv::Mat>> datMA(T, std::vector<cv::Mat>(L));
        datMA = movmean(data);
        int step = std::round(0.5 * opts.cut);
        float maxV = 0;
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        if (data[t][k].at<float>(i,j) > maxV){
                            maxV = data[t][k].at<float>(i,j);
                        }
                    }
                }
            }
        }
        double nSegment = std::max(1.0, std::ceil(T/step)-1);

//        std::cout<<"datMA: "<< std::endl;
//        for (int i = 0; i < 10; ++i) {
//            for (int j = 0; j < 10; ++j) {
//                std::cout<<datMA[0][0].at<float>(i,j)<<" ";
//            }
//            std::cout<<std::endl;
//        }
    }//baselineLinearEstimate()


    void baselineRemoveAndNoiseEstimation(std::vector<std::vector<cv::Mat>>& data1){

        /*
         * smooth the data
         */
//        std::vector<std::vector<cv::Mat>> data1_smo(T);
        std::cout<< "--------start baselineRemoveAndNoiseEstimation--------"<<std::endl;
        int ksize = 2 * ceil(2 * opts.smoXY) + 1;
        std::vector<std::vector<cv::Mat>> data1_smo(T, std::vector<cv::Mat>(L));
        if (opts.smoXY > 0 ){
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    cv::GaussianBlur(data1[t][k],data1_smo[t][k],cv::Size(ksize, ksize), opts.smoXY, opts.smoXY);
                }
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        cv::Mat slice(L,1,CV_32F);
                        for (int k = 0; k < L; ++k) {
                            slice.at<float>(k) = data1_smo[t][k].at<float>(i,j);
                        }
                        cv::GaussianBlur(slice,slice,cv::Size(1,ksize),0,opts.smoXY);
                        for (int k = 0; k < L; ++k) {
                            data1_smo[t][k].at<float>(i,j) = slice.at<float>(k);
                        }
                    }
                }
            }//for(t)
        }//if(smoXY>0)
//        std::cout<<"data1_smo: "<< std::endl;
//        for (int i = 0; i < 10; ++i) {
//            for (int j = 0; j < 10; ++j) {
//                std::cout<<data1_smo[0][0].at<float>(i,j)<<" ";
//            }
//            std::cout<<std::endl;
//        }


        /*
         * linear estimation of F0
         */
        opts.cut = std::min(opts.cut, T);
        //remove baseline
        baselineLinearEstimate(data1_smo);



    }//baselineRemoveAndNoiseEstimation()


}// namespace