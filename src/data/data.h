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
#include <boost/math/distributions/non_central_t.hpp>
#include <omp.h>
#include <fftw3.h>
#include <mat.h>
#include <mex.h>
#include <matrix.h>
#include <cstdlib>
#include <unordered_set>
#include <fstream>

using namespace std;


namespace AQuA{


    struct Point_struct {
        int t;
        int k;
        int i;
        int j;
    };

    struct Score_struct {
        float z_score1;
        float z_score2;
        float t_score1;
        float t_score2;
    };


    bool isDefault();
    void preSettingInit();
    void optsInit();
//    void rawDataSizeInit();
    void Init();
//    void crop(cv::Mat& image, int bdCrop);
    float medianFunc(float* array, int size);
    vector<vector<cv::Mat>> loadData();
    vector<vector<cv::Mat>> load4D_8U(const char* fileName, const char* varName);
    vector<vector<cv::Mat>> load4D(const char* fileName, const char* varName);
    cv::Mat load2D(const char* fileName, const char* varName);
    vector<cv::Mat> load3D(const char* fileName, const char* varName);
    vector<vector<int>> loadCell_int(const char* fileName, const char* varName);
    vector<vector<double>> loadCell_double(const char* fileName, const char* varName);
    vector<cv::Mat> loadCell_matrix(const char* fileName, const char* varName);
    float*** create3dMatrix_float(int h, int w, int l);
    bool*** createEvtSpatialMask(int H, int W, int L);
    int*** create3dMatrix_int(int h, int w, int l);
    int**** create4dMatrix_int(int h, int w, int l, int t);
    void release3dMatrix(float***& data, int h, int w);
    void release3dMatrix_bool(bool***& data, int h, int w);
    void release3dMatrix_int(int***& data, int h, int w);
    mxArray* cvDataToMxArray(const vector<vector<cv::Mat>>& data);
    mxArray* cvDataToMxArray(const vector<cv::Mat>& data);
    mxArray* cvDataToMxArray(const cv::Mat& data);
    void writeDataToMatFile(const vector<vector<cv::Mat>>& data, const string& filename);
    void writeDataToMatFile(const vector<cv::Mat>& data, const string& filename);
    void writeDataToMatFile(cv::Mat& data, const string& filename);
    int sub2ind(int i, int j, int k, int h, int w);
    int sub2ind(int i, int j, int k,int t, int h, int w, int l);
    Point_struct ind2sub(int ind, int h, int w);
    Point_struct ind2sub(int ind, int h, int w, int l);
    void save_vector(const vector<double>& vec, const string& file_path);
    vector<double> load_vector(const string& file_path);





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
         float thrARScl;
         int minDur;
         int minSize;
         int maxSize;
         float circularityThr;
         float spaMergeDist;
         bool needTemp;
         float sigThr;
         float spaSmo;
         float maxDelay;
         bool needRefine;
         bool needGrow;
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
         float step;
         float seedSzRatio;
         int maxSpaScale;
         int minSpaScale;
         int TPatch;
         vector<float> scaleRatios;
         vector<cv::Mat> stdMapOrg1;
         vector<cv::Mat> stdMapGau1;
         vector<cv::Mat> tempVarOrg1;
         vector<cv::Mat> correctPars1;
         vector<vector<cv::Mat>> dF1;
         vector<vector<cv::Mat>> data1_org;
         vector<vector<Point_struct>> arLst1;
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
