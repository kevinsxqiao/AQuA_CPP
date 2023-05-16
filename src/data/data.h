//
// Created by Kevin Qiao on 2/17/2023.
//

#ifndef AQUA_CPP_RAWDATA_H
#define AQUA_CPP_RAWDATA_H


typedef float DATA_TYPE;

#define H AQuA::rawDataSize.size1
#define W AQuA::rawDataSize.size2
#define L AQuA::rawDataSize.size3
#define T AQuA::rawDataSize.frame
#define H_ext (2*H-1)
#define W_ext (2*W-1)
#define L_ext (2*L-1)


#include <string>
#include <iostream>


namespace AQuA{

    bool isDefault();
    void preSettingInit();
    void optsInit();
    void Init();
    void rawDataSizeInit();
    void releaseData(DATA_TYPE**** data, int I, int J, int k);
    void releaseData(double*** data, int I, int J);
    void releaseData(float*** data, int I, int J);
    void releaseData(DATA_TYPE* data);
    void releaseData(int* data);
    DATA_TYPE**** loadData();
    DATA_TYPE**** createSpace();
    DATA_TYPE*** create3dMatrix();
    double*** create3dMatrix_ext_double();
    float*** create3dMatrix_ext_float();
    DATA_TYPE**** create4dMatrix();



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
         DATA_TYPE smoXY;
         DATA_TYPE noiseEstimation;
         int thrARScl;
         int minDur;
         int minSize;
         int maxSize;
         DATA_TYPE circularityThr;
         int needTemp;
         DATA_TYPE sigThr;
         DATA_TYPE maxDelay;
         int needRefine;
         int needGrow;
         int needSpa;
         int sourceSize;
         int cDelay;
         int gtwSmo;
         int ignoreTau;
         DATA_TYPE ratio;
         int regMaskGap;
         int usePG;
         int cut;
         int movAvgWin;
         DATA_TYPE minShow1;
         DATA_TYPE minShowEvtGUI;
         int correctTrend;
         DATA_TYPE propthrmin;
         DATA_TYPE propthrstep;
         DATA_TYPE propthrmax;
         DATA_TYPE compress;
         int gapExt;
         DATA_TYPE frameRate;
         DATA_TYPE spatialRes;
         DATA_TYPE varEst;
         DATA_TYPE fgFluo;
         DATA_TYPE bgFluo;
         int northx;
         int northy;
         int preset;
         std::string filePath1;
         std::string fileName1;
         std::string fileType1;
         std::string filePath2;
         std::string fileName2;
         std::string fileType2;
         double maxValueDat1;
         double minValueDat1;
         int sz[4];
         int BitDepth;
         int maxdF1;
         int maxdF2;
         int singleChannel;
         int alreadyBleachCorrect;
         int alreadyPreprocess;
         int enableTap;
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
