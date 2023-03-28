//
// Created by Kevin Qiao on 3/28/2023.
//

#include "../data/data.h"
#include "regCrossCorrelation.h"
#include <iostream>
namespace AQuA{

    void preProcessRun(){
        if(isDefault() || opts::alreadyProprecess==0 ){ // judge whether this step is already done, since this is time-consuming
            regCrossCorrelation()
        }// if(preSetting_class::isDefault()||...)
    }// preProcessRun()

}

