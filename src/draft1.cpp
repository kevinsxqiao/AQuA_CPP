#include "data/data.h"
#include <iostream>

namespace AQuA{

    void baselineRemoveAndNoiseEstimation(std::vector<std::vector<cv::Mat>>& data1){

        /*
         * smooth the data
         */

    }
}

int main(){
    AQuA::Init();
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> data = AQuA::loadData();
    AQuA::baselineRemoveAndNoiseEstimation(data);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    return 1;
}