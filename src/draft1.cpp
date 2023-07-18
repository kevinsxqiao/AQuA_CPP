#include "data/data.h"
#include "actRun/actRun.h"
#include "preProcessRun/preProcessRun.h"



int main(){
    AQuA::Init();
    auto start = std::chrono::high_resolution_clock::now();
//    std::vector<std::vector<cv::Mat>> dataOrg = AQuA::loadData();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
//    bool*** evtSpatialMask = AQuA::createEvtSpatialMask();
    start = std::chrono::high_resolution_clock::now();
//    AQuA::preProcessRun(dataOrg);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    AQuA::actRun();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
//    AQuA::writeDataToMatFile(dataNew, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test2.mat");


//    AQuA::release3dMatrix_bool(evtSpatialMask,H,W);
    return 0;
}