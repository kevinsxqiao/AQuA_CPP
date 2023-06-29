#include "data/data.h"
#include "preProcess/preProcessRun.h"




int main(){
    AQuA::Init();
    std::vector<std::vector<cv::Mat>> dataOrg = AQuA::loadData();
    auto start = std::chrono::high_resolution_clock::now();
    bool*** evtSpatialMask = AQuA::createEvtSpatialMask();
    std::vector<std::vector<cv::Mat>> dataNew = AQuA::baselineRemoveAndNoiseEstimation(dataOrg, evtSpatialMask);
    AQuA::writeDataToMatFile(dataNew, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test2.mat");
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    AQuA::release3dMatrix_bool(evtSpatialMask,H,W);
    return 0;
}