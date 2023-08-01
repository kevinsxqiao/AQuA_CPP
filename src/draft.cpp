#include "data/data.h"
#include "actRun/actRun.h"
#include "preProcessRun/preProcessRun.h"
#include "data/data.h"
#include "phaseRun/phaseRun.h"

int main(){
    vector<vector<cv::Mat>> mat(5,vector<cv::Mat>(5));
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            mat[i][j] = cv::Mat::zeros(5,5,CV_8U);
        }
    }
//    cv::Mat mat = cv::Mat::zeros(5,5,CV_8U);
//    mat[0][0].at<uchar>(0, 1) = 1;
//    mat[0][0].at<uchar>(2, 2) = 1;
//    mat[4][4].at<uchar>(4, 4) = 1;
    for (int i = 0; i < mat[0][0].rows; ++i) {
        for (int j = 0; j < mat[0][0].cols; ++j) {
            std::cout<<static_cast<int>(mat[0][0].at<uchar>(i,j))<<"  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<endl;
//    cv::Mat mask1 = ((mat == 3) | (mat ==2));
    vector<vector<AQuA::Point_struct>> res = AQuA::bw2Reg(mat);
//    mat.setTo(INFINITY,mask1);
    cout<<res.size()<<" ";
    cout<<res[0].size()<<" ";
//    cv::Mat dst = AQuA::myResize(mat,2,2);
//    for (int i = 0; i < dst.rows; ++i) {
//        for (int j = 0; j < dst.cols; ++j) {
//            std::cout<<dst.at<float>(i,j)<<"  ";
//        }
//        std::cout<<std::endl;
//    }
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


