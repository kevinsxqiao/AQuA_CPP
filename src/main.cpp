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
#include "data/data.h"


namespace AQuA{

    int main() {

        // load the data
        DATA_TYPE**** datOrg1 = loadData();
        DATA_TYPE**** datOrg2 = loadData();
        dataInit();


        int ixTab = menu(); // select action
        switch (ixTab) {
            case 1:
                preProcessRun(datOrg1,datOrg2);
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

        releaseData(datOrg1);
        releaseData(datOrg2);

        return 0;
    }

}
