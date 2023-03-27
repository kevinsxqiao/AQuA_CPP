//
// Created by Kevin Qiao on 2/16/2023.
//

#ifndef AQUA_TEST_H
#define AQUA_TEST_H


namespace AQuA{
    void test111(float (&data1)[rawDataSize::size1][rawDataSize::size2][rawDataSize::size3][rawDataSize::frame]){
        data1[0][0][0][0] = 0;
    }
    void testC(float (&c)[2 * rawDataSize::size1 - 1][2 * rawDataSize::size2 - 1][2 * rawDataSize::size3 - 1]){
        for(int i=0;i<AQuA::rawDataSize::size1;++i){
            for (int j = 0; j < AQuA::rawDataSize::size2; ++j) {
                for (int k = 0; k < AQuA::rawDataSize::size3; ++k) {
                    c[i][j][k] = 1.5;
                }
            }
        }
    }
}

#endif //AQUA_TEST_H
