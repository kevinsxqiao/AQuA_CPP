//
// Created by Kevin Qiao on 3/28/2023.
//

#include "preProcessRun.h"

namespace AQuA{

    void preProcessRun(DATA_TYPE**** data1, DATA_TYPE**** data2){
        if(isDefault() || opts::alreadyPreprocess==0 ){ // judge whether this step is already done, since this is time-consuming
            DATA_TYPE**** datOrg1 = regCrossCorrelation(data1,data2);
            /*
             *     scl.hrg = [1,size(datOrg1,1)];
                    scl.wrg = [1,size(datOrg1,2)];
                    scl.lrg = [1,size(datOrg1,3)];
             */

        }
    }// preProcessRun()

}

