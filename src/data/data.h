//
// Created by Kevin Qiao on 2/17/2023.
//

#ifndef AQUA_CPP_RAWDATA_H
#define AQUA_CPP_RAWDATA_H

#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>
#include <vector>
#include <boost/math/distributions/normal.hpp>
#include <omp.h>
#include <fftw3.h>
#include <mat.h>
#include <mex.h>
#include <matrix.h>


#define H AQuA::rawDataSize.size1
#define W AQuA::rawDataSize.size2
#define L AQuA::rawDataSize.size3
#define T AQuA::rawDataSize.frame
#define H_ext (2*H-1)
#define W_ext (2*W-1)
#define L_ext (2*L-1)


namespace AQuA{

    bool isDefault();
    void preSettingInit();
    void optsInit();
    void rawDataSizeInit();
    void Init();
//    void crop(cv::Mat& image, int bdCrop);
    std::vector<std::vector<cv::Mat>> loadData();
    float*** create3dMatrix(int h, int w, int l);
    void release3dMatrix(float***& data, int h, int w);
    void release3dMatrix_bool(bool***& data, int h, int w);
    mxArray* cvDataToMxArray(const std::vector<std::vector<cv::Mat>>& data);
    void writeDataToMatFile(std::vector<std::vector<cv::Mat>>& data, const std::string& filename);
//    float**** createSpace();
//    float*** create3dMatrix();
//    double*** create3dMatrix_ext_double();
//    float*** create3dMatrix_ext_float();
//    float**** create4dMatrix();
//    void releaseData(float**** data, int I, int J, int k);
//    void releaseData(double*** data, int I, int J);
//    void releaseData(float*** data, int I, int J);
//    void releaseData(float* data);
//    void releaseData(int* data);



    struct rawDataSize_struct{ // clarification structure of raw data size
        int size1;
        int size2;
        int size3;
        int frame;
    };// struct

    extern rawDataSize_struct rawDataSize;



    struct preSetting_struct{//  member clarification
         const int registrateCorrect_default = 0;
         const int bleachCorrect_default = 0;
         int registrateCorrect;
         int bleachCorrect;
    };//struct

    extern preSetting_struct preSetting;


    struct opts_struct{
         int registrateCorrect;
         int bleachCorrect;
         float smoXY;
         float noiseEstimation;
         int thrARScl;
         int minDur;
         int minSize;
         int maxSize;
         float circularityThr;
         int needTemp;
         float sigThr;
         float maxDelay;
         int needRefine;
         int needGrow;
         int needSpa;
         int sourceSize;
         int cDelay;
         int gtwSmo;
         int ignoreTau;
         float ratio;
         int regMaskGap;
         int usePG;
         int cut;
         int movAvgWin;
         float minShow1;
         float minShowEvtGUI;
         int correctTrend;
         float propthrmin;
         float propthrstep;
         float propthrmax;
         float compress;
         int gapExt;
         float frameRate;
         float spatialRes;
         float varEst;
         float fgFluo;
         float bgFluo;
         int northx;
         int northy;
         int preset;
         const char* fileName1;
         const char* fileName2;
         double maxValueDat1;
         double minValueDat1;
         int sz[4];
         int BitDepth;
         float maxdF1;
         float maxdF2;
         int singleChannel;
         bool alreadyBleachCorrect;
         bool alreadyPreprocess;
         int enableTap;
         double medSmo;
         std::vector<cv::Mat> stdMapOrg1;
         std::vector<cv::Mat> stdMapGau1;
         std::vector<cv::Mat> tempVarOrg1;
         std::vector<cv::Mat> correctPars1;
    };// struct opts

    extern opts_struct opts;


    struct scl_struct{
         int min;
         int max;
         int bri1;
         int bri2;
         int briL;
         int briR;
         float briOv;
         float minOv;
         float maxOv;
         int hrg[2];
         int wrg[2];
         int lrg[2];
         int h;
         int w;
         int l;
    };// struct scl

    extern scl_struct scl;


}// namespace
#endif //AQUA_CPP_RAWDATA_H
