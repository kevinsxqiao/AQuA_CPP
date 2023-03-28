//
// Created by Kevin Qiao on 2/10/2023.
//

#include <iostream>
#include "preProcess/preProcessRun.h"
#include "actRun.h"
#include "phaseRun.h"
#include "evtRun.h"
#include "feaRun.h"
#include "flow/flow.h"



int main() {
    // load raw data


    int ixTab = menu(); // select action
    switch (ixTab) {
        case 1:
            AQuA::preProcessRun();
            break;
//        case 2:
//            AQuA::actRun();
//            break;
//        case 3:
//            AQuA::phaseRun();
//            break;
//        case 4:
//            AQuA::evtRun();
//            break;
//        case 5:
//            AQuA::feaRun();
//            break;
    }
    return 0;
}