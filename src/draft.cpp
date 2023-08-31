#include "data/data.h"
#include "actRun/actRun.h"
#include "preProcessRun/preProcessRun.h"
#include "data/data.h"
#include "phaseRun/phaseRun.h"


int main() {

    vector<cv::Mat> covMatrixs = AQuA::loadCell_matrix("../cfg/Order_mus_sigmas.mat","covMatrixs");
    cout<<covMatrixs[0].size()<<endl;
    cout<<covMatrixs[1].size()<<endl;
    cout<<covMatrixs[2].size()<<endl;
    cout<<covMatrixs[1].at<float>(0,0)<<endl;
    cout<<covMatrixs[1].at<float>(0,1)<<endl;

    return 0;
}





