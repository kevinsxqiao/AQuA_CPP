#include "data/data.h"
#include <iostream>
//#include <opencv2/opencv.hpp>
//#include "preProcess/baselineRemoveAndNoiseEstimation.h"

namespace AQuA{




}// namespace

int main(){
    AQuA::Init();
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> data = AQuA::loadData();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    return 1;
}