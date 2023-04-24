#include <opencv2/opencv.hpp>
#include "data/data.h"

int main() {

    std::vector<cv::Mat> images;
    std::string img_name = "C:/Users/Kevin Qiao/Desktop/3D_data/3D_dataFrame ";
    std::string name_exd = ".tif";

    extern AQuA::rawDataSize_struct rawDataSizeStruct;
    AQuA::rawDataSizeInit();
    int frameNumber = 100;
    for (int t = 1; t < frameNumber; ++t) {
        cv::imreadmulti(img_name + std::to_string(t) + name_exd, images, cv::IMREAD_GRAYSCALE);
        if(t==1){
            AQuA::rawDataSize.size1 = images[0].rows;
            W = images[0].cols;
            H = images[0].rows;
            L = images.size();
        }
        for (int i = 0; i < L; ++i) {
            cv::imshow("length", images[i]);
            cv::waitKey();
        }
    }



    return 0;
}
