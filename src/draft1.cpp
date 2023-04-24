#include <opencv2/opencv.hpp>
#include "data/data.h"

void loadData() {

    std::vector<cv::Mat> images;
    std::string pre_name = "C:/Users/Kevin Qiao/Desktop/3D_data/3D_dataFrame ";
    std::string name_exd = ".tif";

    AQuA::rawDataSizeInit();
    T = 100;
    AQuA::create4dMatrix();

    for (int t = 1; t < T; ++t) {
        cv::imreadmulti(pre_name + std::to_string(t) + name_exd, images, cv::IMREAD_GRAYSCALE);
        if (t == 1) {
            W = images[0].cols;
            H = images[0].rows;
            L = images.size();
        }
        for (int i = 0; i < L; ++i) {
            cv::imshow("length", images[i]);
            std::cout << "length " << i + 1 << " of " << L << " frame " << t << " of " << T << std::endl;
            cv::waitKey();
        }
    }
}

int main(){
    loadData();
    return 0;
}

