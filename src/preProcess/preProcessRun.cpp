//
// Created by Kevin Qiao on 3/28/2023.
//

#include "preProcessRun.h"

namespace AQuA{

    void preProcessRun(DATA_TYPE**** data1){
        if(isDefault() || opts::alreadyProprecess==0 ){ // judge whether this step is already done, since this is time-consuming
            regCrossCorrelation();
        }
    }// preProcessRun()

}

