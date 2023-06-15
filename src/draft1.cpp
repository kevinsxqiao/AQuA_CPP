#include "data/data.h"
#include <iostream>
//#include <cmath>
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


    std::vector<cv::Mat> truncated_kept_var(const std::vector<cv::Mat>& quantiles){
        std::vector<cv::Mat> pars(L);
        boost::math::normal_distribution<float> normal_dist;
        for (int k = 0; k < L; ++k) {
            pars[k] = cv::Mat(H,W,CV_32F);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    float quantile = quantiles[k].at<float>(i,j);
                    if (quantile == 0.0f){
                        pars[k].at<float>(i,j) = 2.0f;
                    }
                    else{
                        float a = boost::math::quantile(normal_dist, quantile);
                        float phi_a = boost::math::pdf(normal_dist,a);
                        float mu = a * quantile + phi_a;
                        float second_order = a * a * quantile + 1 - quantile + a * phi_a;
                        pars[k].at<float>(i,j) = (second_order - mu * mu) * 2;
                    }
                }
            }
        }//for(k)

        return pars;
    }//truncated_kept_var()


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
//        std::cout<<"datMA: "<< std::endl;
////        for (int k = 0; k < 2; ++k) {
//            for (int i = 0; i < 7; ++i) {
//                for (int j = 0; j < 7; ++j) {
//                    std::cout<<datMA[0][0].at<float>(i,j)<<" ";
//                }
//                std::cout<<std::endl;
//            }
//        }

        int step = static_cast<int>(std::round(0.5 * opts.cut));
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
//        std::vector<std::vector<cv::Mat>> minPosition(nSegment, std::vector<cv::Mat>(L, cv::Mat::ones(H, W, CV_32S)));
        std::vector<std::vector<cv::Mat>> minPosition(nSegment, std::vector<cv::Mat>(L));
        for (int kk = 0; kk < nSegment; ++kk) {
            int t0 = kk * step;
            int t1 = std::min(T, t0 + opts.cut);
            for (int k = 0; k < L; ++k) {
                minPosition[kk][k] = cv::Mat::ones(H, W, CV_32S);
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        float minV = 10;
                        for (int t = t0; t < t1; ++t) {
                            if (datMA[t][k].at<float>(i, j) < minV) {
                                minV = datMA[t][k].at<float>(i, j);
                                minPosition[kk][k].at<int>(i, j) = t;
                            }//if
                        }//for(t)
                    }
                }
            }//for(k)
//            float minV = 10;
//            for (int t = t0; t < t1; ++t) {
//                if (datMA[t][0].at<float>(1,1) < minV) {
//                    minV = datMA[t][0].at<float>(1,1);
//                    minPosition[kk][0].at<int>(1,1) = t;
//                }
//            }
        }//for(kk)

//        for (int k = 0; k < 10; ++k) {
//            std::cout<<k<<" :"<< minPosition[0][k].at<int>(1,1)<<"  "<< std::endl;
//        }
//        std::cout<<"minPosition"<<std::endl;
//        std::cout<<minPosition[0][0].at<int>(1,1)<<" : ";
//        std::cout<<datMA[minPosition[0][0].at<int>(1,1)][0].at<float>(1,1)<<"   ";
//        std::cout<<std::endl;
//        std::cout<<"datMA[t]"<<std::endl;
//        for (int t = 0; t < 200; ++t) {
//            std::cout<<t<<" :"<< datMA[t][0].at<float>(1,2)<<"  ";
//        }
//        std::cout<<"minPosition"<<std::endl;
//        for (int i = 1; i < 20; ++i) {
//            std::cout<<minPosition[0][0].at<int>(i,1)<<"  ";
//        }
//        std::cout<<std::endl;


        std::vector<std::vector<cv::Mat>> F0(T, std::vector<cv::Mat>(L));
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                F0[t][k] = cv::Mat::ones(H, W, CV_32F);
            }
        }
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
                            value.push_back(datMA[curP[tt]][k].at<float>(i,j));
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
                            F0[t][k].at<float>(i,j) = curve.at<float>(0,t);
                        }
                    }//for(j)
                }//for(i)
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

//        std::cout<<"F0: "<< std::endl;
//        std::cout<<"height of image:"<< F0[0][0].rows << std::endl;
//        std::cout<<"width of image:"<< F0[0][0].cols << std::endl;
//        std::cout<<"length of image:"<< F0[0].size() << std::endl;
//        std::cout<<"time frames of image:"<< F0.size() << std::endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                std::cout<<F0[0][0].at<float>(i,j)<<" ";
//            }
//            std::cout<<std::endl;
//        }

        return F0;
    }//baselineLinearEstimate()


    void noiseEstimationFunction(std::vector<std::vector<cv::Mat>>& dataOrg){
        /*
         * variance map
         * calculate the variance of raw data
         */
        bool correctNoise = true;
        std::vector<cv::Mat> tempMap(L);
        std::vector<cv::Mat> tempVarOrg(L);
        for (int k = 0; k < L; ++k) {
            tempMap[k] = cv::Mat(H,W,CV_32F);
            tempVarOrg[k] = cv::Mat(H,W,CV_32F);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    double sum = 0;
                    for (int t = 0; t < T-1; ++t) {
                        sum += pow(dataOrg[t][k].at<float>(i,j) - dataOrg[t+1][k].at<float>(i,j), 2);
                    }
                    tempMap[k].at<float>(i,j) = static_cast<float>(sum/T);
                }
            }
            tempVarOrg[k] = tempMap[k] / 2;
        }//for(k)
        if (correctNoise){
            std::vector<cv::Mat> countInValid(L);
            std::vector<cv::Mat> totalSamples(L);
            std::vector<cv::Mat> ratio(L);
            for (int k = 0; k < L; ++k) {
                tempMap[k] = cv::Mat(H,W,CV_32F);
                tempVarOrg[k] = cv::Mat(H,W,CV_32F);
                ratio[k] = cv::Mat(H,W,CV_32F);
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        int count_invalid = 0, count_samples=0;
                        for (int t = 0; t < T; ++t) {
                            if (dataOrg[t][k].at<float>(i,j) == 0){
                                ++count_invalid;
                            }
                            if (!isnan(dataOrg[t][k].at<float>(i,j))){
                                ++count_samples;
                            }
                            countInValid[k].at<int>(i,j) = count_invalid;
                            totalSamples[k].at<int>(i,j) = count_samples;
                        }
                    }
                }
                ratio[k] = countInValid[k] / totalSamples[k];
            }//for(k)
            std::vector<cv::Mat> correctPars = truncated_kept_var(ratio);
            std::vector<cv::Mat> varMapOrg(L);
            for (int k = 0; k < L; ++k) {
                varMapOrg[k] = cv::Mat(H,W,CV_32F);
                varMapOrg[k] = tempMap[k] / correctPars[k];
            }
        }//if(correctNoise)
        else{
            std::vector<cv::Mat> varMapOrg(L);
            for (int k = 0; k < L; ++k) {
                varMapOrg[k] = cv::Mat(H,W,CV_32F);
                varMapOrg[k] = tempVarOrg[k];
            }
        }//else

    }//noiseEstimationFunction()


    void baselineRemoveAndNoiseEstimation(std::vector<std::vector<cv::Mat>>& dataOrg){
        /*
         * smooth the data
         */
        std::cout<< "--------start baselineRemoveAndNoiseEstimation--------"<<std::endl;
        int ksize = 2 * ceil(2 * opts.smoXY) + 1;
        std::vector<std::vector<cv::Mat>> dataSmo(T, std::vector<cv::Mat>(L));
        if (opts.smoXY > 0 ){
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    cv::GaussianBlur(dataOrg[t][k],dataSmo[t][k],cv::Size(ksize, ksize), opts.smoXY, opts.smoXY);
                }
//                // 3rd dimension smoothing
//                for (int i = 0; i < H; ++i) {
//                    for (int j = 0; j < W; ++j) {
//                        cv::Mat slice(L,1,CV_32F);
//                        for (int k = 0; k < L; ++k) {
//                            slice.at<float>(k) = dataSmo[t][k].at<float>(i,j);
//                        }
//                        cv::GaussianBlur(slice,slice,cv::Size(1,ksize),0,opts.smoXY);
//                        for (int k = 0; k < L; ++k) {
//                            dataSmo[t][k].at<float>(i,j) = slice.at<float>(k);
//                        }
//                    }
//                }
            }//for(t)
        }//if(smoXY>0)

//        std::cout<<"dataSmo: "<< std::endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                std::cout<<dataSmo[0][0].at<float>(i,j)<<" ";
//            }
//            std::cout<<std::endl;
//        }


        /*
         * linear estimation of F0
         */
        opts.cut = std::min(opts.cut, T);
        //remove baseline
        std::vector<std::vector<cv::Mat>> F0(T, std::vector<cv::Mat>(L));
        std::vector<cv::Mat> F0Pro(L);
        for (int k = 0; k < L; ++k) {
            F0Pro[k] = cv::Mat(H,W,CV_32F);
        }
        F0 = baselineLinearEstimate(dataSmo);
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
//        noiseEstimationFunction();



    }//baselineRemoveAndNoiseEstimation()


}// namespace

int main(){
    auto start = std::chrono::high_resolution_clock::now();
    AQuA::Init();
    std::vector<std::vector<cv::Mat>> data = AQuA::loadData();
    AQuA::baselineRemoveAndNoiseEstimation(data);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    return 1;
}