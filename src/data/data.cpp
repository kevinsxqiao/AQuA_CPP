//
// Created by Kevin Qiao on 3/28/2023.
//

#include "data.h"


namespace AQuA{


    DATA_TYPE**** create4dMatrix(){

        DATA_TYPE ****data;
        data = new DATA_TYPE ***[H];
        for (int i = 0; i < H; ++i) {
            data[i] = new DATA_TYPE**[W];

        }
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                data[i][j] = new DATA_TYPE*[L];
            }
        }
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    data[i][j][k] = new DATA_TYPE[T];
                }
            }
        }

        int value = 1;

//        assign values
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


    void releaseData(DATA_TYPE**** data){
        delete data;
    }// releaseData()


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