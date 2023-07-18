//
// Created by Kevin Qiao on 7/18/2023.
//

#include "phaseRun.h"

namespace AQuA{

    void phaseRun(){
        std::cout<< "--------start temporal segmentation--------"<<std::endl;
        if (opts.needTemp){
            opts.step = 0.5; // sigma
            opts.sigThr = 3.5; // Z score of seed significance
            opts.maxDelay = 0.6; // allowed maximum dissimilarity in merging
            opts.seedSzRatio = 0.01; // seed size relative to active region
            opts.needRefine = false; // peaks are temporally adjacent, need regine
            opts.needGrow = false; // grow active regions according to signal pattern


        }

    }

}
