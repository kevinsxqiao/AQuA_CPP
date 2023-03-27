//
// Created by Kevin Qiao on 2/10/2023.
//

#ifndef AQUA_PREPROCESSRUN_CPP
#define AQUA_PREPROCESSRUN_CPP
#include "data/data.h"


namespace AQuA{
    void preProcessRun(){ // activate voxels detection and update overlay map
        if(preSetting_class::isDefault() || opts_struct::alreadyProprecess==0){ // judge whether this step is already done, since this is time-consuming
            std::cout << 1;
        }// if(preSetting_class::isDefault()||...)
    }// preProcessRun()
}// namespace

#endif //AQUA_PREPROCESSRUN_CPP
