#include "data/data.h"
#include "preProcess/regCrossCorrelation.h"


int main(){
    AQuA::Init();
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> data = AQuA::loadData();
//    AQuA::baselineRemoveAndNoiseEstimation(data);
    AQuA::regCrossCorrelation(data);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    return 1;
}