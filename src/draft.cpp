#include "data/data.h"
#include "actRun/actRun.h"
#include "preProcessRun/preProcessRun.h"
#include "data/data.h"
#include "phaseRun/phaseRun.h"

int main(){
    cv::Mat mat = cv::Mat::ones(5,5,CV_32F);
    mat.at<float>(0, 1) = 2.0f;
    mat.at<float>(1, 1) = 3.0f;
    mat.at<float>(2, 2) = 3.0f;
    mat.at<float>(4, 4) = 3.0f;
    mat.at<float>(4, 0) = 3.0f;
    for (int i = 0; i < mat.rows; ++i) {
        for (int j = 0; j < mat.cols; ++j) {
            std::cout<<mat.at<float>(i,j)<<"  ";
        }
        std::cout<<std::endl;
    }
    cv::Mat dst = AQuA::myResize(mat,2,2);
    for (int i = 0; i < dst.rows; ++i) {
        for (int j = 0; j < dst.cols; ++j) {
            std::cout<<dst.at<float>(i,j)<<"  ";
        }
        std::cout<<std::endl;
    }
//    std::cout<<mat.type()<<std::endl;
//    mat.at<uchar>(1,1) = 0;
//    mat.at<uchar>(2,2) = 0;
//    int flag = 0;
//    for (int i = 0; i < 5; ++i) {
//        for (int j = 0; j < 5; ++j) {
//            std::cout<<static_cast<int>(mat.at<uchar>(i,j))<<"  ";
//            if (mat.at<uchar>(i,j) !=1){
//                flag++;
//            }
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<<flag<<std::endl;
//    std::cout<<"size of 8U: "<<sizeof (uchar)<<std::endl;
//    std::cout<<"size of 8S: "<<sizeof (char)<<std::endl;
//    std::cout<<"size of 32S: "<<sizeof (int)<<std::endl;

    return 0;
}


