#include "data/data.h"
#include "actRun/actRun.h"



int main(){
    std::vector<std::vector<cv::Mat>> mat(5,std::vector<cv::Mat>(5));
    for (int t = 0; t < 5; ++t) {
        for (int k = 0; k < 5; ++k) {
            mat[t][k] = cv::Mat::zeros(5,5,CV_32S);
        }
    }

    mat[0][0].at<int>(0,0) = 1;
//    mat[0][0].at<int>(1,0) = 1;
//    mat[0][1].at<int>(1,2) = 1;

    mat[0][1].at<int>(4,4) = 1;
    mat[0][1].at<int>(3,4) = 1;
//    mat[0][2].at<int>(1,2) = 1;
//    mat[0][2].at<int>(1,0) = 1;

    mat[4][0].at<int>(1,1) = 1;
    std::vector<std::vector<AQuA::Point_struct>> poi = AQuA::bw2Reg(mat);
    std::cout<< poi.size()<<std::endl;
    std::cout<< poi[0].size()<<std::endl;
    std::cout<< poi[1].size()<<std::endl;
    std::cout<< poi[1][0].t<< poi[1][0].k<< poi[1][0].i<< poi[1][0].j<<std::endl;


    return 0;
}