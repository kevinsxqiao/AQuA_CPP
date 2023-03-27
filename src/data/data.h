//
// Created by Kevin Qiao on 2/17/2023.
//

#ifndef AQUA_CPP_RAWDATA_H
#define AQUA_CPP_RAWDATA_H

namespace AQuA{


    struct rawDataSize{ // clarification structure of raw data size
        static const short size1;
        static const short size2;
        static const short size3;
        static const short frame;
    };// struct

     const short rawDataSize::size1 = 2;// static member definition
     const short rawDataSize::size2 = 2;
     const short rawDataSize::size3 = 2;
     const short rawDataSize::frame = 2;

    float datOrg1[rawDataSize::size1][rawDataSize::size2][rawDataSize::size3][rawDataSize::frame] = {1,2,3,4,5,
                                                                                                     6,7,8,9,10,
                                                                                                     11,12,13,14,15,16};
    float datOrg2[rawDataSize::size1][rawDataSize::size2][rawDataSize::size3][rawDataSize::frame] = {0};


    struct preSetting_structure{// static member clarification
        static const short registrateCorrect_default = 0;
        static const short bleachCorrect_default = 0;
        static short registrateCorrect;
        static short bleachCorrect;
    };//struct




    // preSetting_structure static member definition
    short preSetting_structure::registrateCorrect = preSetting_structure::registrateCorrect_default;
    short preSetting_structure::bleachCorrect = preSetting_structure::bleachCorrect_default;





    struct opts_struct{
        static short registrateCorrect;
        static short bleachCorrect;
        static float smoXY;
        static float noiseEstimation;
        static short thrARScl;
        static short minDur;
        static int minSize;
        static int maxSize;
        static float circularityThr;
        static short needTemp;
        static float sigThr;
        static float maxDelay;
        static short needRefine;
        static short needGrow;
        static short needSpa;
        static int sourceSize;
        static short cDelay;
        static short gtwSmo;
        static short ignoreTau;
        static float ratio;
        static short regMaskGap;
        static short usePG;
        static short cut;
        static short movAvgWin;
        static float minShow1;
        static float minShowEvtGUI;
        static short correctTrend;
        static float propthrmin;
        static float propthrstep;
        static float propthrmax;
        static float compress;
        static short gapExt;
        static float frameRate;
        static float spatialRes;
        static float varEst;
        static float fgFluo;
        static float bgFluo;
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
        static short alreadyProprecess;
        static short enableTap;
    };// struct opts


    // opts_struct static member definition
    short opts_struct::alreadyProprecess = 0;


    bool isDefault() { //judge if the value is changed or not; ---- true = no change; false = changed
        return (preSetting_structure::registrateCorrect == preSetting_structure::registrateCorrect_default
                && preSetting_structure::bleachCorrect == preSetting_structure::bleachCorrect_default);
    }


}// namespace
#endif //AQUA_CPP_RAWDATA_H
