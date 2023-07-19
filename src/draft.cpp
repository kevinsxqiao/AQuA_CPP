#include "data/data.h"
#include "actRun/actRun.h"
#include "preProcessRun/preProcessRun.h"
#include "data/data.h"

int main(){
//    std::vector<std::vector<cv::Mat>> mat(5,std::vector<cv::Mat>(5));
    cv::Mat temp1 = cv::Mat::zeros(5,5,CV_32S);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> mat2(5000);
    for (int t = 0; t < 5000; ++t) {
        for (int k = 0; k < 5000; ++k) {
            mat2[t].emplace_back(temp1.clone());
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> mat1(5000);
    for (int t = 0; t < 5000; ++t) {
        mat1[t].resize(5000);
        for (int k = 0; k < 5000; ++k) {
            mat1[t][k] = mat2[t][k];
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> mat3(5000);
    for (int t = 0; t < 5000; ++t) {
        mat3[t].resize(5000);
        for (int k = 0; k < 5000; ++k) {
            mat3[t][k] = mat2[t][k].clone();
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    std::cout<<&mat1<<" "<<&mat2<<std::endl;
    std::cout<<&mat3<<" "<<&mat2<<std::endl;
//    start = std::chrono::high_resolution_clock::now();
//    std::vector<std::vector<cv::Mat>> mat(5000,std::vector<cv::Mat>(5000));
//    for (int t = 0; t < 5000; ++t) {
//        for (int k = 0; k < 5000; ++k) {
//            mat[t][k] = temp1.clone();
//        }
//    }
//    end = std::chrono::high_resolution_clock::now();
//    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
//    cv::Mat temp;
//    cv::multiply(mat[0][0],2,temp);
//    for (int i = 0; i < 5; ++i) {
//        for (int j = 0; j < 5; ++j) {
//            std::cout<<mat[3][2].at<int>(i,j)<<" ";
//        }
//        std::cout<<std::endl;
//    }
//std::cout<<mat.size()<<std::endl;
//std::cout<<mat[0].size()<<std::endl;
//
//    mat[0][0].at<int>(0,0) = 1;
//    mat[1][0].at<int>(0,0) = 1;
//
//    mat[1][0].at<int>(1,2) = 1;
//
//    mat[0][1].at<int>(4,4) = 1;
//    mat[0][1].at<int>(3,4) = 1;
////    mat[0][2].at<int>(1,2) = 1;
////    mat[0][2].at<int>(1,0) = 1;
//
//    mat[4][0].at<int>(1,1) = 1;
////    std::cout<<AQuA::sub2ind(0,1,1,0,5,5,5)<<std::endl;
//
//    std::vector<std::vector<AQuA::Point_struct>> poi = AQuA::bw2Reg(mat);
//    for (int i = 0; i < poi.size(); ++i) {
//        for (int j = 0; j < poi[i].size(); ++j) {
//            std::cout<<poi[i][j].t<<" "<<poi[i][j].k<<" "<<poi[i][j].i<<" "<<poi[i][j].j<<" ";
//            std::cout<<"    ";
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<< poi[1][0].t<< poi[1][0].k<< poi[1][0].i<< poi[1][0].j<<std::endl;


    return 0;
}


