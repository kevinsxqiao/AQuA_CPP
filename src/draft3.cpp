#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <mat.h>
#include <vector>
#include <boost/math/distributions/normal.hpp>
#include "data/data.h"
#include <omp.h>



int main(){
    omp_set_num_threads(12);
    auto start = std::chrono::high_resolution_clock::now();
    AQuA::Init();
    std::vector<std::vector<cv::Mat>> dataOrg = AQuA::loadData();
    std::cout<< "0,0"<<std::endl;
    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            std::cout<< dataOrg[0][0].at<float>(i,j)<<"  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<< "0,1"<<std::endl;
    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            std::cout<< dataOrg[0][1].at<float>(i,j)<<"  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<< "1,1"<<std::endl;
    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            std::cout<< dataOrg[1][1].at<float>(i,j)<<"  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    return 0;
}