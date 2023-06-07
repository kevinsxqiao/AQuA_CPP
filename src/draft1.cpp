#include "data/data.h"
#include <iostream>

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


    void baseLineLinearEstimate(std::vector<std::vector<cv::Mat>>& data){
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

//        std::cout<<"datMA: "<< std::endl;
//        for (int i = 0; i < 10; ++i) {
//            for (int j = 0; j < 10; ++j) {
//                std::cout<<datMA[0][0].at<float>(i,j)<<" ";
//            }
//            std::cout<<std::endl;
//        }
    }//baseLineLinearEstimate()


    void baselineRemoveAndNoiseEstimation(std::vector<std::vector<cv::Mat>>& data1){

        /*
         * smooth the data
         */
//        std::vector<std::vector<cv::Mat>> data1_smo(T);
        int ksize = 2 * ceil(2 * opts.smoXY) + 1;
        std::vector<std::vector<cv::Mat>> data1_smo(T, std::vector<cv::Mat>(L));
//        for (int t = 0; t < T; ++t) {
//            for (int k = 0; k < L; ++k) {
//                data1_smo[t].emplace_back(data1[t][k].clone());
//            }
//        }
        if (opts.smoXY > 0 ){
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    cv::GaussianBlur(data1[t][k],data1_smo[t][k],cv::Size(ksize, ksize), opts.smoXY, opts.smoXY);
                }
            }
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
        baseLineLinearEstimate(data1_smo);



    }//baselineRemoveAndNoiseEstimation()

}

int main(){
    AQuA::Init();
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> data = AQuA::loadData();
    AQuA::baselineRemoveAndNoiseEstimation(data);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    return 1;
}