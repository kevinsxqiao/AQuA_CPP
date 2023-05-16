#include <opencv2/opencv.hpp>
#include "data/data.h"


void crop(cv::Mat& image, int bdCrop){
//    cv::Rect roi = cv::Rect(bdCrop,bdCrop,W-2*bdCrop,H-2*bdCrop);
//    image.adjustROI(-roi.y, image.rows - (roi.y+roi.height), -roi.x, image.cols- (roi.x+roi.width));
    image.adjustROI(-bdCrop, -bdCrop, -bdCrop, -bdCrop);
//    std::cout<<"height of image:"<< image.rows << std::endl;
//    std::cout<<"width of image:"<< image.cols << std::endl;
}

DATA_TYPE **** loadData() {

    T = 100;
    std::vector<cv::Mat> frame[T];
    std::string pre_name = "C:/Users/Kevin Qiao/Desktop/3D_data/3D_dataFrame ";
    std::string name_ext = ".tif";
    DATA_TYPE min,max=0;
    DATA_TYPE mmin=255, mmax=0;
    int bdCrop = AQuA::opts.regMaskGap;
    int BitDepth = -1;


//    AQuA::rawDataSizeInit();

//    DATA_TYPE **** data = AQuA::create4dMatrix();
    std::cout<<"reading data ..."<<std::endl;
    for (int t = 0; t < 1; ++t) {
        cv::imreadmulti(pre_name + std::to_string(t+1) + name_ext, frame[t],  cv::IMREAD_GRAYSCALE); // uchar 8-bit unsigned int
        if (t == 0) {
            H = frame[0][0].rows;
            W = frame[0][0].cols;
            L = static_cast<int>(frame[0].size());
            std::cout<<"height of image:"<< H << std::endl;
            std::cout<<"width of image:"<< W << std::endl;
            std::cout<<"length of image:"<< L << std::endl;
        }//(t=0) output size of each dimension

        for (int k = 0; k < L; ++k) {
//            for (int i = 0; i < H; ++i) {
//                for (int j = 0; j < W; ++j) {
//////                    frame[t][k].at<uchar>(i,j) = 255 -  frame[t][k].at<uchar>(i,j);
////                    std::cout<<static_cast<int>(frame[t][k].at<uchar>(i,j))<<std::endl;
//                }//for(j)
//            }//for(i)
            crop(frame[t][k],bdCrop);
            cv::minMaxLoc(frame[t][k], reinterpret_cast<double *>(&min), reinterpret_cast<double *>(&max));
            if (min < mmin){
                mmin = min;
            }
            if(max > mmax){
                mmax = max;
            }
//            std::cout<< min<< " "<<max<< std::endl;
//            cv::imwritemulti("C:/Users/Kevin Qiao/Desktop/1.tif",frame[t]); //save image
            std::cout << "slice " << k + 1 << " of " << L << " frame " << t + 1 << " of " << T << std::endl;
            cv::imshow("length", frame[t][k]);
            cv::waitKey();
        }//for(k)
    }//for(t)

    AQuA::opts.maxValueDat1 = mmax;
    AQuA::opts.minValueDat1 = mmin;

//    for (int t = 0; t < T; ++t) {
//        for (int k = 0; k < L; ++k) {
//            for (int i = 0; i < H; ++i) {
//                for (int j = 0; j < W; ++j) {
//                    frame[t][k].at<uchar>(i,j) -= AQuA::opts.minValueDat1;
//                    if (frame[t][k].at<uchar>(i,j)<0){
//                        frame[t][k].at<uchar>(i,j) = 0;
//                    }
//                    frame[t][k].at<uchar>(i,j) = frame[t][k].at<uchar>(i,j) / (AQuA::opts.maxValueDat1 - AQuA::opts.minValueDat1);
//                }//for(i)
//            }//for(j)
//            std::cout << "slice " << k + 1 << " of " << L << " frame " << t + 1 << " of " << T << std::endl;
//            cv::imshow("length", frame[t][k]);
//            cv::waitKey();
//        }//for(k)
//    }//for(t)

    AQuA::opts.sz[0] = H;
    AQuA::opts.sz[1] = W;
    AQuA::opts.sz[2] = L;
    AQuA::opts.sz[3] = T;
    AQuA::opts.BitDepth = BitDepth;

}//loadData()



int main(){
    AQuA::rawDataSizeInit();
    AQuA::optsInit();
    AQuA::preSettingInit();
    loadData();
    return 0;
}

