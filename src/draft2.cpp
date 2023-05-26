#include <mat.h>
#include <iostream>
#include <string>
#include <vector>
#include <opencv2/opencv.hpp>
#include "data/data.h"

namespace AQuA{
    /*
     *  read the data from mat file and store in std::vector<std::vector<cv::Mat>> frame(T)
     *  access pixel value by frame[t][k].at<float>(i,j)
     */
    std::vector<std::vector<cv::Mat>> loadData_1() {
        MATFile *pmatFile;
        mxArray *pMxArray;
        double *pdata;
        int bdCrop = AQuA::opts.regMaskGap;
        double min,max=0;
        double mmin=255, mmax=0;
        int BitDepth = -1;
        float normalizedParameter;

        std::cout<< "--------loading data--------"<<std::endl;
        const char *filename = "C:/Users/Kevin Qiao/Desktop/AQuA_data/Test_global_local_3D.mat";
        pmatFile = matOpen(filename, "r");
        if (pmatFile == nullptr) {
            std::cout<< "--------error opening file--------"<<std::endl;
            std::exit(-1);
        }

        pMxArray = matGetVariable(pmatFile, "synthetic");
        if (pMxArray == nullptr) {
            std::cout<< "--------error reading variable from file--------"<<std::endl;
            std::exit(-1);
        }

        pdata = mxGetPr(pMxArray);
        if (pdata == nullptr) {
            std::cout<< "--------error reading data from variable-------"<<std::endl;
            std::exit(-1);
        }

        const mwSize *dims = mxGetDimensions(pMxArray);
        H = dims[0];
        W = dims[1];
        L = dims[2];
        T = dims[3];
//        std::cout<<"original size: "<< std::endl;
//        std::cout<<"height of image:"<< H << std::endl;
//        std::cout<<"width of image:"<< W << std::endl;
//        std::cout<<"length of image:"<< L << std::endl;
//        std::cout<<"time frames of image:"<< T << std::endl;

        std::vector<std::vector<cv::Mat>> frame(T);
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                frame[t].emplace_back(H- 2*bdCrop,W- 2*bdCrop,CV_32F);
                for (int i = 0, i_src = bdCrop; i < H - 2*bdCrop; ++i, ++ i_src) {
                    for (int j = 0, j_src= bdCrop; j < W - 2*bdCrop; ++j, ++j_src) {
                        frame[t][k].at<float>(i,j) = static_cast<float>(pdata[j_src*H + i_src + k*H*W + t*H*W*L]);
                    }//for(j)
                }//for(i)
                cv::minMaxLoc(frame[t][k], &min, &max);
                if (min < mmin){
                    mmin = min;
                }
                if(max > mmax){
                    mmax = max;
                }
            }//for(k)
        }//for(t)

        H = dims[0] - 2*bdCrop;
        W = dims[1] - 2*bdCrop;
        std::cout<<"after cropping: "<< std::endl;
        std::cout<<"height of image:"<< H << std::endl;
        std::cout<<"width of image:"<< W << std::endl;
        std::cout<<"length of image:"<< L << std::endl;
        std::cout<<"time frames of image:"<< T << std::endl;

        //release MAT pointer
        if (pMxArray != nullptr) {
            mxDestroyArray(pMxArray);
        }

        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }

        AQuA::opts.maxValueDat1 = mmax;
        AQuA::opts.minValueDat1 = mmin;
//        std::cout<<"minValue: "<< opts.minValueDat1 << "  maxValue: "<<opts.maxValueDat1<<std::endl;
        normalizedParameter = static_cast<float>(mmax -mmin);
        std::cout<<"--------data loaded--------"<<std::endl;

        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        frame[t][k].at<float>(i,j) = (frame[t][k].at<float>(i,j) - static_cast<float>(AQuA::opts.minValueDat1)) / normalizedParameter;
//                    std::cout<< frame[t][k].at<float>(i,j)<< "  ";// display pixel value
                    }//for(i)
                }//for(j)
            }//for(k)
        }//for(t)
        AQuA::opts.sz[0] = H;
        AQuA::opts.sz[1] = W;
        AQuA::opts.sz[2] = L;
        AQuA::opts.sz[3] = T;
        AQuA::opts.BitDepth = BitDepth;

//        for (int i = 0; i < 10; ++i) {
//            for (int j = 0; j < 10; ++j) {
//                std::cout<< frame[2][3].at<float>(i,j)<< " ";
//            }
//            std::cout<<std::endl;
//        }

        return frame;
    }//loadData()
}//namespace


int main(){
    AQuA::Init();
    auto start = std::chrono::high_resolution_clock::now();
    AQuA::loadData_1();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    return 1;
}
