//
// Created by Kevin Qiao on 2/10/2023.
//

#include "preProcessRun/preProcessRun.h"
#include "actRun/actRun.h"
#include "phaseRun.h"
#include "evtRun.h"
#include "feaRun.h"
#include "flow/flow.h"
#include "data/data.h"

    int main() {

        AQuA::Init();
        std::vector<std::vector<cv::Mat>> data1 = AQuA::loadData();
        int ixTab = AQuA::menu(); // select action
        while(1){
            switch (ixTab) {
                case 0:
                    return 0;
                case 1:
                    AQuA::preProcessRun(data1);
                    break;
                case 2:
                    AQuA::actRun();
                    break;
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
        }

    }


