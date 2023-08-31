#include "data/data.h"
//#include "actRun/actRun.h"
//#include "preProcessRun/preProcessRun.h"
#include "phaseRun/phaseRun.h"

using namespace std;

namespace AQuA {
}






int main(){
    std::cout << std::fixed << std::setprecision(4);
    AQuA::Init();
    auto start = std::chrono::high_resolution_clock::now();
//    std::vector<std::vector<cv::Mat>> dataOrg = AQuA::loadData();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
//    bool*** evtSpatialMask = AQuA::createEvtSpatialMask();
    start = std::chrono::high_resolution_clock::now();
//    std::vector<std::vector<cv::Mat>> dF1 = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/Test_global_local_3D.mat","synthetic");
//        std::vector<std::vector<int>> arlst = AQuA::loadCell("C:/Users/Kevin Qiao/Desktop/AQuA_data/arlst.mat","arLst1");
//    AQuA::seDetection(dF1,datOrg1,arLst1);
    AQuA::phaseRun();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
//    AQuA::writeDataToMatFile(dataNew, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test2.mat");


//    AQuA::release3dMatrix_bool(evtSpatialMask,H,W);
    return 0;
}