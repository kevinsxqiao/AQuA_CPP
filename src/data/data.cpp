//
// Created by Kevin Qiao on 3/28/2023.
//

#include "data.h"


namespace AQuA{

    /*
     * define some data structure
     */
    struct rawDataSize_struct rawDataSize;
    struct preSetting_struct preSetting;
    struct opts_struct opts;

    /*
     * create a 4d matrix, which size is T,L,H,W
     * initialize the values with 0
     */
    DATA_TYPE**** create4dMatrix(){

        DATA_TYPE**** data;
        data = new DATA_TYPE*** [H];
        for (int i = 0; i < H; ++i) {
            data[i] = new DATA_TYPE** [W];

        }
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                data[i][j] = new DATA_TYPE* [L];
            }
        }
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    data[i][j][k] = new DATA_TYPE [T];
                }
            }
        }

//        for (int i = 0; i < H; ++i) {
//            for (int j = 0; j < W; ++j) {
//                for (int k = 0; k < L; ++k) {
//                    for (int t = 0; t < T; ++t) {
//                        data[i][j][k][t] = 0;
//                    }
//                }
//            }
//        }


        return data;
    }// getData()


    /*
     * create a 3d matrix, which size is H,W,L
     * initialize the values with 0
     */
    DATA_TYPE*** create3dMatrix(){
        DATA_TYPE*** data;
        data = new DATA_TYPE** [H];
        for (int i = 0; i < H; ++i) {
            data[i] = new DATA_TYPE* [W];
        }
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                data[i][j] = new DATA_TYPE [L];
            }
        }
//        for (int i = 0; i < H; ++i) {
//            for (int j = 0; j < W; ++j) {
//                for (int k = 0; k < L; ++k) {
//                    data[i][j][k] = 0;
//                }
//            }
//        }
        return data;
    }// create3dMatrix()


    /*
     * create a 3d matrix which size is H_ext,W_ext,L_ext
     * initialize the values with 0
     */
    double *** create3dMatrix_ext_double(){
        double*** data;
        data = new double ** [H_ext];
        for (int i = 0; i < H_ext; ++i) {
            data[i] = new double * [W_ext];
        }
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                data[i][j] = new double [L_ext];
            }
        }
//        for (int i = 0; i < H_ext; ++i) {
//            for (int j = 0; j < W_ext; ++j) {
//                for (int k = 0; k < L_ext; ++k) {
//                    data[i][j][k] = 0;
//                }
//            }
//        }
        return data;
    }// create3dMatrix_ext()


    float *** create3dMatrix_ext_float(){
        float *** data;
        data = new float ** [H_ext];
        for (int i = 0; i < H_ext; ++i) {
            data[i] = new float * [W_ext];
        }
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                data[i][j] = new float [L_ext];
            }
        }
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                for (int k = 0; k < L_ext; ++k) {
                    data[i][j][k] = 0;
                }
            }
        }
        return data;
    }// create3dMatrix_ext()_float()


//    load the file and return a 4-d matrix pointer
    DATA_TYPE**** loadData(){
        /*
         * load tiff file
         */
        H = 2;
        W = 2;
        L = 2;
        T = 2;
        return create4dMatrix();
    }// loadData()


    void releaseData(int* data){
        delete[] data;
        data = NULL;
    }// releaseData

    void releaseData(DATA_TYPE* data){
        delete[] data;
        data = NULL;
    }// releaseData

    void releaseData(DATA_TYPE**** data, int I, int J, int K){
        for (int i = 0; i < I; ++i) {
            for (int j = 0; j < J; ++j) {
                for (int k = 0; k < K; ++k) {
                    delete[] data[i][j][k];
                }
                delete[] data[i][j];
            }
            delete[] data[i];
        }
        delete[] data;
        data = NULL;
    }// releaseData


    void releaseData(double*** data, int I, int J){
        for (int i = 0; i < I; ++i) {
            for (int j = 0; j < J; ++j) {
                delete[] data[i][j];
            }
            delete[] data[i];
        }
        delete[] data;
        data = NULL;
    }// releaseData

    void releaseData(float*** data, int I, int J){
        for (int i = 0; i < I; ++i) {
            for (int j = 0; j < J; ++j) {
                delete[] data[i][j];
            }
            delete[] data[i];
        }
        delete[] data;
        data = NULL;
    }// releaseData


//    judge if registration and bleach have been executed; ---- true = no ; false = both executed
    bool isDefault() {
        return (preSetting.registrateCorrect == preSetting.registrateCorrect_default
                && preSetting.bleachCorrect == preSetting.bleachCorrect_default);
    }// isDefault()

    void preSettingInit(){
        preSetting.registrateCorrect = preSetting.registrateCorrect_default;
        preSetting.bleachCorrect = preSetting.bleachCorrect_default;
        std::cout<< "    preSetting initialized! "<<std::endl;
    }// preSettingInit()


    void optsInit(){
        opts.alreadyPreprocess = 0;
        opts.alreadyBleachCorrect = 0;
        opts.movAvgWin = 25;
        opts.cut = 200;
        opts.smoXY = 0.5;
        opts.BitDepth;
        opts.regMaskGap = 5;
        opts.singleChannel = true;
        opts.registrateCorrect = 1;
        opts.bleachCorrect = 1;
        opts.smoXY = 0.5;

        std::cout<< "    opts initialized! "<<std::endl;
    }// optsInit()

    void rawDataSizeInit(){
        rawDataSize={0};
        std::cout<< "    rawDataSize initialized! "<<std::endl;
    }

    void Init(){
        preSettingInit();
        optsInit();
        rawDataSizeInit();
        std::cout<<std::endl;
    }







}// namespace