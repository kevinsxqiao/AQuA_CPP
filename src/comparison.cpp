#include <iostream>
#include <opencv2/opencv.hpp>
#include <chrono>

using namespace cv;

const int NUM_ITERATIONS = 100;

int main() {
    Mat image = imread("C:/Users/Kevin Qiao/Desktop/3D_data/3D_dataFrame 1.tif", IMREAD_GRAYSCALE);
    if (image.empty()) {
        std::cerr << "Failed to read image" << std::endl;
        return -1;
    }

    std::cout<< "execute threshold: "<<std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        Mat processed = image.clone();
        threshold(processed, processed, 128, 255, THRESH_BINARY);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_mat = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using cv::Mat: " << duration_mat << " milliseconds" << std::endl;


    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        Mat processed_1d[2];
        processed_1d[1]= image.clone();
        threshold(processed_1d[1], processed_1d[1], 128, 255, THRESH_BINARY);
    }
    end = std::chrono::high_resolution_clock::now();
    duration_mat = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using cv::Mat[]: " << duration_mat << " milliseconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        Mat processed_2d[2][2];
        processed_2d[1][1]= image.clone();
        threshold(processed_2d[1][1], processed_2d[1][1], 128, 255, THRESH_BINARY);
    }
    end = std::chrono::high_resolution_clock::now();
    duration_mat = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using cv::Mat[][]: " << duration_mat << " milliseconds" << std::endl;


    int rows = image.rows;
    int cols = image.cols;
    unsigned char* data = image.data;
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        unsigned char* processed = new unsigned char[rows * cols];
        for (int j = 0; j < rows * cols; ++j) {
            processed[j] = data[j] > 128 ? 255 : 0;
        }
        delete[] processed;
    }
    end = std::chrono::high_resolution_clock::now();
    auto duration_array = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using Array: " << duration_array << " milliseconds" << std::endl;

    rows = image.rows;
    cols = image.cols;
    unsigned char** data_2d = new unsigned char*[rows];
    for (int i = 0; i < rows; ++i) {
        data_2d[i] = new unsigned char[cols];
        for (int j = 0; j < cols; ++j) {
            data_2d[i][j] = image.at<unsigned char>(i, j);
        }
    }


    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        unsigned char** processed = new unsigned char*[rows];
        for (int j = 0; j < rows; ++j) {
            processed[j] = new unsigned char[cols];
            for (int k = 0; k < cols; ++k) {
                processed[j][k] = data_2d[j][k] > 128 ? 255 : 0;
            }
        }

        for (int j = 0; j < rows; ++j) {
            delete[] processed[j];
        }
        delete[] processed;
    }
    end = std::chrono::high_resolution_clock::now();
    duration_array = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using 2D Array: " << duration_array << " milliseconds" << std::endl;


    for (int i = 0; i < rows; ++i) {
        delete[] data_2d[i];
    }
    delete[] data_2d;
    std::cout<< std::endl;
    std::cout<< "execute pixel+1: "<<std::endl;
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        Mat processed = image.clone();
        for (int j = 0; j < processed.rows; ++j) {
            for (int k = 0; k < processed.cols; ++k) {
                processed.at<unsigned char>(j, k) += 1;
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration_mat = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using cv::Mat: " << duration_mat << " milliseconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        Mat processed_1d1[2];
        processed_1d1[1]= image.clone();
        for (int j = 0; j < processed_1d1[1].rows; ++j) {
            for (int k = 0; k < processed_1d1[1].cols; ++k) {
                processed_1d1[1].at<unsigned char>(j, k) += 1;
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration_mat = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using cv::Mat[]: " << duration_mat << " milliseconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        Mat processed_2d2[2][2];
        processed_2d2[1][1]= image.clone();
        for (int j = 0; j < processed_2d2[1][1].rows; ++j) {
            for (int k = 0; k < processed_2d2[1][1].cols; ++k) {
                processed_2d2[1][1].at<unsigned char>(j, k) += 1;
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration_mat = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using cv::Mat[][]: " << duration_mat << " milliseconds" << std::endl;


    rows = image.rows;
    cols = image.cols;
    unsigned char* data_p = image.data;
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        unsigned char* processed = new unsigned char[rows * cols];
        for (int j = 0; j < rows * cols; ++j) {
            processed[j] = data_p[j] + 1;
        }
        delete[] processed;
    }
    end = std::chrono::high_resolution_clock::now();
    duration_array = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using Array: " << duration_array << " milliseconds" << std::endl;

    rows = image.rows;
    cols = image.cols;

    unsigned char** data_array = new unsigned char* [rows];
    for (int i = 0; i < rows; ++i) {
        data_array[i] = new unsigned char[cols];
        for (int j = 0; j < cols; ++j) {
            data_array[i][j] = image.at<unsigned char>(i, j);
        }
    }

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        unsigned char** processed_array = new unsigned char* [rows];
        for (int j = 0; j < rows; ++j) {
            processed_array[j] = new unsigned char[cols];
            for (int k = 0; k < cols; ++k) {
                processed_array[j][k] = data_array[j][k] + 1;
            }
        }
        for (int j = 0; j < rows; ++j) {
            delete[] processed_array[j];
        }
        delete[] processed_array;
    }
    end = std::chrono::high_resolution_clock::now();
    duration_array = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Using 2D Array: " << duration_array << " milliseconds" << std::endl;


    for (int i = 0; i < rows; ++i) {
        delete[] data_array[i];
    }
    delete[] data_array;

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> mat1(5000);
    cv::Mat temp1 = cv::Mat::zeros(5,5,CV_32S);
    for (int t = 0; t < 5000; ++t) {
        mat1[t].resize(5000);
        for (int k = 0; k < 5000; ++k) {
            mat1[t][k] = temp1.clone();
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> mat2(5000);
    for (int t = 0; t < 5000; ++t) {
        for (int k = 0; k < 5000; ++k) {
            mat2[t].emplace_back(temp1.clone());
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> mat(5000,std::vector<cv::Mat>(5000));
    for (int t = 0; t < 5000; ++t) {
        for (int k = 0; k < 5000; ++k) {
            mat[t][k] = temp1.clone();
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;

    return 0;
}

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