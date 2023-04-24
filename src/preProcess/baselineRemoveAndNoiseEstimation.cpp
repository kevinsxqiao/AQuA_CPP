//
// Created by Kevin Qiao on 3/28/2023.
//

#include "baselineRemoveAndNoiseEstimation.h"


namespace AQuA{
    DATA_TYPE**** baselineRemoveAndNoiseEstimation(DATA_TYPE**** data){

//        % smooth the data (memory-efficient version)
        DATA_TYPE **** datSmo = create4dMatrix();
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    for (int t = 0; t < T; ++t) {
                        datSmo[i][j][k][t] = data[i][j][k][t];
                    }
                }
            }
        }
//3d gaussian filter
        if (opts.smoXY > 0){
            for (int t = 0; t < T; ++t) {

            }
        }

//linear estimation of F0
        opts.cut = std::min(opts.cut, T);
//        remove baseline



        return data;
    }
}// namespace