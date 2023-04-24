#include <opencv2/opencv.hpp>
#include "data/data.h"

DATA_TYPE **** loadData() {

    std::vector<cv::Mat> images;
    std::string pre_name = "C:/Users/Kevin Qiao/Desktop/3D_data/3D_dataFrame ";
    std::string name_exd = ".tif";

//    AQuA::rawDataSizeInit();
    T = 100;
    DATA_TYPE **** data = AQuA::create4dMatrix();

    for (int t = 1; t < T; ++t) {
        cv::imreadmulti(pre_name + std::to_string(t) + name_exd, images, cv::IMREAD_GRAYSCALE);
        if (t == 1) {
            H = images[0].rows;
            W = images[0].cols;
            L = images.size();
        }

        for (int k = 0; k < L; ++k) {
            std::cout<< H << "  "<< W << "  "<< L<< std::endl;
            std::cout<< images[k].type()<< "  ";
            std::cout<<images[k].depth()<< "  ";
            std::cout<<images[k].channels()<< "  ";
            images[k].convertTo(images[k],CV_32F, 1.0/255.0);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
//                    data[i][j][k][t] = images[k].at<float>(i,j); //uchar:: unsigned 8-bit integer
//                    std::cout<<data[i][j][k][t]<<std::endl;
                    std::cout<<images[k].at<float>(i,j)<<std::endl;
                }//for(j)
            }//for(i)

            cv::imshow("length", images[k]);

            std::cout<< images[k].type()<< "  ";
            std::cout<<images[k].depth()<< "  ";
            std::cout<<images[k].channels()<< "  ";
            std::cout << "length " << k + 1 << " of " << L << " frame " << t << " of " << T << std::endl;
            cv::waitKey();
        }//for(k)
    }//for(t)
}

int main(){
    loadData();
    return 0;
}

