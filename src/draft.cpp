#include "data/data.h"
#include <iostream>
#include <opencv2/opencv.hpp>

namespace AQuA{

}

int main(){
        cv::Mat mat  = cv::Mat::zeros(7,7,CV_32F);;
        mat.at<float>(3,3) = 1;
        for (int i = 0; i < 7; ++i) {
            for (int j = 0; j < 7; ++j) {
                std::cout<<mat.at<float>(i,j)<< "  ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
        cv::GaussianBlur(mat,mat,cv::Size(5,5),1,1);
        for (int i = 0; i < 7; ++i) {
            for (int j = 0; j < 7; ++j) {
                std::cout<< std::fixed << std::setprecision(4)<<mat.at<float>(i,j)<< "  ";
            }
            std::cout<<std::endl;
        }
}