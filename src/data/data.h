//
// Created by Kevin Qiao on 2/17/2023.
//

#ifndef AQUA_CPP_RAWDATA_H
#define AQUA_CPP_RAWDATA_H

#define DATA_TYPE double

#define H AQuA::rawDataSize::size1
#define W AQuA::rawDataSize::size2
#define L AQuA::rawDataSize::size3
#define T AQuA::rawDataSize::frame
#define H_ext (2*H-1)
#define W_ext (2*W-1)
#define L_ext (2*L-1)


#include <string>
#include <iostream>


namespace AQuA{


    struct rawDataSize{ // clarification structure of raw data size
        static short size1;
        static short size2;
        static short size3;
        static short frame;
    };// struct



    struct preSetting{// static member clarification
        static const short registrateCorrect_default = 0;
        static const short bleachCorrect_default = 0;
        static short registrateCorrect;
        static short bleachCorrect;
    };//struct


    struct opts{
        static short registrateCorrect;
        static short bleachCorrect;
        static DATA_TYPE smoXY;
        static DATA_TYPE noiseEstimation;
        static short thrARScl;
        static short minDur;
        static int minSize;
        static int maxSize;
        static DATA_TYPE circularityThr;
        static short needTemp;
        static DATA_TYPE sigThr;
        static DATA_TYPE maxDelay;
        static short needRefine;
        static short needGrow;
        static short needSpa;
        static int sourceSize;
        static short cDelay;
        static short gtwSmo;
        static short ignoreTau;
        static DATA_TYPE ratio;
        static short regMaskGap;
        static short usePG;
        static short cut;
        static short movAvgWin;
        static DATA_TYPE minShow1;
        static DATA_TYPE minShowEvtGUI;
        static short correctTrend;
        static DATA_TYPE propthrmin;
        static DATA_TYPE propthrstep;
        static DATA_TYPE propthrmax;
        static DATA_TYPE compress;
        static short gapExt;
        static DATA_TYPE frameRate;
        static DATA_TYPE spatialRes;
        static DATA_TYPE varEst;
        static DATA_TYPE fgFluo;
        static DATA_TYPE bgFluo;
        static short northx;
        static short northy;
        static short preset;
        static std::string filePath1;
        static std::string fileName1;
        static std::string fileType1;
        static std::string filePath2;
        static std::string fileName2;
        static std::string fileType2;
        static int maxValueDat1;
        static int minValueDat1;
        static short sz[3];
        static short BitDepth;
        static short maxdF1;
        static short maxdF2;
        static short singleChannel;
        static short alreadyBleachCorrect;
        static short alreadyPreprocess;
        static short enableTap;
    };// struct opts

    struct scl{
        static short min;
        static short max;
        static short bri1;
        static short bri2;
        static short briL;
        static short briR;
        static float briOv;
        static float minOv;
        static float maxOv;
        static short hrg[2];
        static short wrg[2];
        static short lrg[2];
        static short h;
        static short w;
        static short l;
    };// struct scl



    bool isDefault();
    void dataInit();
    void releaseData(DATA_TYPE**** data);
    DATA_TYPE**** loadData();
    DATA_TYPE**** createSpace();


}// namespace
#endif //AQUA_CPP_RAWDATA_H
