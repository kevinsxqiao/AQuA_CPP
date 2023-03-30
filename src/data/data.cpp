//
// Created by Kevin Qiao on 3/28/2023.
//

#include "data.h"


namespace AQuA{

    /*
     * create a 4d matrix, which size is H,W,L,T
     * initialize the values with 1++
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


        int value = 1;
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    for (int t = 0; t < T; ++t) {
                        data[i][j][k][t] = value++;
                    }
                }
            }
        }



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
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    data[i][j][k] = 0;
                }
            }
        }
        return data;
    }// create3dMatrix()


    /*
     * create a 3d matrix which size is H_ext,W_ext,L_ext
     * initialize the values with 0
     */
    DATA_TYPE*** create3dMatrix_ext(){
        DATA_TYPE*** data;
        data = new DATA_TYPE** [H_ext];
        for (int i = 0; i < H_ext; ++i) {
            data[i] = new DATA_TYPE* [W_ext];
        }
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                data[i][j] = new DATA_TYPE [L_ext];
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
    }// create3dMatrix_ext()


    float*** create3dMatrix_ext_float(){
        float*** data;
        data = new float** [H_ext];
        for (int i = 0; i < H_ext; ++i) {
            data[i] = new float* [W_ext];
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
        T  = 2;
        return create4dMatrix();
    }// loadData()


    void releaseData(int* data){
        delete data;
        data = NULL;
    }// releaseData

    void releaseData(DATA_TYPE**** data){
        delete data;
        data = NULL;
    }// releaseData

    void releaseData(DATA_TYPE* data){
        delete data;
        data = NULL;
    }// releaseData

    void releaseData(DATA_TYPE*** data){
        delete data;
        data = NULL;
    }// releaseData

    void releaseData(float*** data){
        delete data;
        data = NULL;
    }// releaseData


//    judge if registration and bleach have been executed; ---- true = no ; false = both executed
    bool isDefault() {
        return (preSetting::registrateCorrect == preSetting::registrateCorrect_default
                && preSetting::bleachCorrect == preSetting::bleachCorrect_default);
    }// isDefault()

    void dataInit(){
        preSetting::registrateCorrect = preSetting::registrateCorrect_default;
        preSetting::bleachCorrect = preSetting::bleachCorrect_default;
        opts::alreadyPreprocess = 0;
        opts::alreadyBleachCorrect = 0;
        std::cout<< " preSetting and opts initialized! "<<std::endl;
    }


}// namespace