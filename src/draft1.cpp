#include "data/data.h"
#include <iostream>
//#include <opencv2/opencv.hpp>
//#include "preProcess/baselineRemoveAndNoiseEstimation.h"

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


    std::vector<std::vector<cv::Mat>> baselineLinearEstimate(std::vector<std::vector<cv::Mat>>& data){
        std::vector<std::vector<cv::Mat>> datMA(T, std::vector<cv::Mat>(L));
        datMA = movmean(data);
//        for (int t = 0; t < T; t++) {
//            for (int k = 0; k < L; k++) {
//                cv::blur(data[t][k], datMA[t][k], cv::Size(opts.movAvgWin, opts.movAvgWin));
//            }
//        }
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        if (isnan(data[t][k].at<float>(i,j))){
                            datMA[t][k].at<float>(i,j) = NAN;
                        }
                    }
                }
            }
        }
        std::cout<<"datMA: "<< std::endl;
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                std::cout<<datMA[0][0].at<float>(i,j)<<" ";
            }
            std::cout<<std::endl;
        }

        int step = std::round(0.5 * opts.cut);
        float maxV = 0;
        double maxVal = 0;
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                cv::minMaxLoc(cv::Mat(data[t][k]), nullptr, &maxVal);
                if (maxVal > maxV){
                    maxV = static_cast<float>(maxVal);
                }
            }
        }
        int nSegment = static_cast<int>(std::max(1.0, std::ceil(T/step)-1));
        std::vector<std::vector<cv::Mat>> minPosition(nSegment, std::vector<cv::Mat>(L, cv::Mat::ones(H, W, CV_32F)));
        for (int kk = 0; kk < nSegment; ++kk) {
            int t0 = kk * step;
            int t1 = std::min(T, t0 + opts.cut);
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        float minV = 1.0;
                        for (int t = t0; t < t1; ++t) {
                            if (datMA[t][k].at<float>(i, j) < minV) {
                                minV = datMA[t][k].at<float>(i, j);
                                minPosition[kk][k].at<int>(i, j) = t;
                            }
                        }
                    }
                }
            }//for(k)
        }
        for (int i = 0; i < 10; ++i) {
            std::cout<<minPosition[0][0].at<int>(0,i)<<" "<<std::endl;
        }
            std::vector<std::vector<cv::Mat>> F0(T, std::vector<cv::Mat>(L, cv::Mat::zeros(H, W, CV_32F)));
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        std::vector<int> curP;
                        std::vector<float> value;
                        for (int kk = 0; kk < nSegment; ++kk) {
                            curP.push_back(minPosition[kk][k].at<int>(i,j));
                        }
                        if (nSegment > 1){
                            std::sort(curP.begin(),curP.end());
                            auto it = std::unique(curP.begin(),curP.end());
                            curP.erase(it, curP.end());
                        }
                        for (int tt = 0; tt < curP.size(); ++tt) {
                            value.push_back(datMA[tt][k].at<float>(i,j));
//                            if (!isnan(value[tt])){
//                                curP[tt] = curP[tt];
//                                value[tt] = value[tt];
//                            }
                        }

                        int nMin = value.size();
                        cv::Mat curve = cv::Mat::zeros(1,T,CV_32F);
                        if (nMin == 0){
                            curve = maxV;
                        }
                        else{
                            //first part
                            for (int l = 0; l < curP[0]; ++l) {
                                curve.at<float>(0,l) = value[0];
                            }
                            //end part
                            for (int l = curP[nMin-1]; l < T; ++l) {
                                curve.at<float>(0,l) = value[nMin-1];
                            }
                            //middle part
                            for (int l = 0; l < nMin-1; ++l) {
                                int mt1 = curP[l];
                                int mt2 = curP[l+1];
                                for (int m = mt1, index=0; m < mt2; ++m, ++index) {
                                    curve.at<float>(0,m) = value[l] + (value[l+1]-value[l])/static_cast<float>((mt2-mt1)*index);
                                }

                            }
                        }//else
                        for (int t = 0; t < T; ++t) {
                            F0[t][k].at<float>(i,j) = curve.at<float>(t);
                        }
                    }//for(j)
                }
            }//for(k)

        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        if (isnan(data[t][k].at<float>(i,j))){
                            F0[t][k].at<float>(i,j) = maxV;
                        }
                    }
                }
            }
        }

        std::cout<<"F0: "<< std::endl;
        std::cout<<"height of image:"<< F0[0][0].rows << std::endl;
        std::cout<<"width of image:"<< F0[0][0].cols << std::endl;
        std::cout<<"length of image:"<< F0[0].size() << std::endl;
        std::cout<<"time frames of image:"<< F0.size() << std::endl;
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                std::cout<<F0[0][0].at<float>(i,j)<<" ";
            }
            std::cout<<std::endl;
        }

        return F0;
    }//baselineLinearEstimate()


    void baselineRemoveAndNoiseEstimation(std::vector<std::vector<cv::Mat>>& data1){
        /*
         * smooth the data
         */
        std::cout<< "--------start baselineRemoveAndNoiseEstimation--------"<<std::endl;
        int ksize = 2 * ceil(2 * opts.smoXY) + 1;
        std::vector<std::vector<cv::Mat>> data1_smo(T, std::vector<cv::Mat>(L));
        if (opts.smoXY > 0 ){
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    cv::GaussianBlur(data1[t][k],data1_smo[t][k],cv::Size(ksize, ksize), opts.smoXY, opts.smoXY);
                }
//                // 3rd dimension smoothing
//                for (int i = 0; i < H; ++i) {
//                    for (int j = 0; j < W; ++j) {
//                        cv::Mat slice(L,1,CV_32F);
//                        for (int k = 0; k < L; ++k) {
//                            slice.at<float>(k) = data1_smo[t][k].at<float>(i,j);
//                        }
//                        cv::GaussianBlur(slice,slice,cv::Size(1,ksize),0,opts.smoXY);
//                        for (int k = 0; k < L; ++k) {
//                            data1_smo[t][k].at<float>(i,j) = slice.at<float>(k);
//                        }
//                    }
//                }
            }//for(t)
        }//if(smoXY>0)

        std::cout<<"data1_smo: "<< std::endl;
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                std::cout<<data1_smo[0][0].at<float>(i,j)<<" ";
            }
            std::cout<<std::endl;
        }


        /*
         * linear estimation of F0
         */
        opts.cut = std::min(opts.cut, T);
        //remove baseline
        std::vector<std::vector<cv::Mat>> F0(T, std::vector<cv::Mat>(L));
        std::vector<cv::Mat> F0Pro(L,cv::Mat(H,W,CV_32F));
        F0 = baselineLinearEstimate(data1_smo);
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    float sum = 0;
                    for (int t = 0; t < T; ++t) {
                        sum += F0[t][k].at<float>(i,j);
                    }
                    F0Pro[k].at<float>(i,j) = sum/T;
                }
            }
        }

        //noise estimation




    }//baselineRemoveAndNoiseEstimation()


}// namespace

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