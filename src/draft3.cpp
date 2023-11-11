#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <mat.h>
#include <vector>
#include <boost/math/distributions/normal.hpp>
#include "data/data.h"
#include <omp.h>
#include "phaseRun/phaseRun.h"
#include <chrono>




#include <opencv2/opencv.hpp>
#include <vector>
#include <iostream>
#include <chrono>

int main() {
    int dim1 = 50;  // 第一维大小
    int dim2 = 200;  // 第二维大小
    int totalMats = dim1 * dim2;
    int x = 50, y = 50; // 特定位置
    float sum;
    std::chrono::steady_clock::time_point start, end;


    // 使用 vector<vector<cv::Mat>>
    std::vector<std::vector<cv::Mat>> matsVec(dim1, std::vector<cv::Mat>(dim2, cv::Mat::ones(100, 100, CV_32F)));

    start = std::chrono::steady_clock::now();
    sum = 0;
    for (const auto& row : matsVec) {
        for (const auto& mat : row) {
            sum += mat.at<float>(y, x);
        }
    }
    end = std::chrono::steady_clock::now();
    std::cout << "Vector time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;


    // 使用一维数组
    cv::Mat* mats1D = new cv::Mat[totalMats];
    for (int i = 0; i < totalMats; ++i) {
        mats1D[i] = cv::Mat::ones(100, 100, CV_32F);
    }

    start = std::chrono::steady_clock::now();
    sum = 0;
    for (int i = 0; i < totalMats; ++i) {
        sum += mats1D[i].at<float>(y, x);
    }
    end = std::chrono::steady_clock::now();
    std::cout << "One-dimensional array time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;
//    delete[] mats1D;

    // 使用二维数组
    cv::Mat** mats2D = new cv::Mat*[dim1];
    for (int i = 0; i < dim1; ++i) {
        mats2D[i] = new cv::Mat[dim2];
        for (int j = 0; j < dim2; ++j) {
            mats2D[i][j] = cv::Mat::ones(100, 100, CV_32F);
        }
    }

    start = std::chrono::steady_clock::now();
    sum = 0;
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            sum += mats2D[i][j].at<float>(y, x);
        }
    }
    end = std::chrono::steady_clock::now();
    std::cout << "Two-dimensional array time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;
//    for (int i = 0; i < dim1; ++i) {
//        delete[] mats2D[i];
//    }
//    delete[] mats2D;



    std::vector<std::vector<cv::Mat>> matsVec1(dim1,vector<cv::Mat>(dim2));
    for (int i = 0; i < dim1; ++i) {
        for (int j = 0; j < dim2; ++j) {
            matsVec1[i][j] = cv::Mat::ones(100,100,CV_32F);
        }
    }

    start = std::chrono::steady_clock::now();
    sum = 0;
    for (const auto& row : matsVec1) {
        for (const auto& mat : row) {
            sum += mat.at<float>(y, x);
        }
    }
    end = std::chrono::steady_clock::now();
    std::cout << "Vector1 time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;

    return 0;
}
