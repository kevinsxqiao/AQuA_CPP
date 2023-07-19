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
