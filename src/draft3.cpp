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
    double min,max=0;
    double mmin=255, mmax=0;
    int bdCrop = AQuA::opts.regMaskGap;
    int BitDepth = -1;
    DATA_TYPE normalizedParameter;


//    DATA_TYPE **** data = AQuA::create4dMatrix();
    std::cout<<"reading data ..."<<std::endl;
    for (int t = 0; t < T; ++t) {
        cv::imreadmulti(pre_name + std::to_string(t+1) + name_ext, frame[t],  cv::IMREAD_GRAYSCALE); // uchar 8-bit unsigned int
        if (t == 0) {
            H = frame[0][0].rows;
            W = frame[0][0].cols;
            L = static_cast<int>(frame[0].size());
//            std::cout<<"after cropping:"<< std::endl;
//            std::cout<<"height of image:"<< H << std::endl;
//            std::cout<<"width of image:"<< W << std::endl;
//            std::cout<<"length of image:"<< L << std::endl;
        }//(t=0) output size of each dimension
        for (int k = 0; k < L; ++k) {
//            for (int i = 0; i < H; ++i) {
//                for (int j = 0; j < W; ++j) {
//////                    frame[t][k].at<uchar>(i,j) = 255 -  frame[t][k].at<uchar>(i,j);
////                    std::cout<<static_cast<int>(frame[t][k].at<uchar>(i,j))<<std::endl;
//                }//for(j)
//            }//for(i)
            crop(frame[t][k],bdCrop);
            if (k == 0) {
                H = frame[0][0].rows;
                W = frame[0][0].cols;
            std::cout<<"  after cropping:"<< std::endl;
            std::cout<<"height of image:"<< H << std::endl;
            std::cout<<"width of image:"<< W << std::endl;
            std::cout<<"length of image:"<< L << std::endl;
            }//(t=0) output size of each dimension
            cv::minMaxLoc(frame[t][k], &min, &max);
            if (min < mmin){
                mmin = min;
            }
            if(max > mmax){
                mmax = max;
            }
//            std::cout<< min<< " "<<max<< std::endl;
//            cv::imwritemulti("C:/Users/Kevin Qiao/Desktop/1.tif",frame[t]); //save image
//            std::cout << "slice " << k + 1 << " of " << L << " frame " << t + 1 << " of " << T << std::endl;
//            cv::imshow("length", frame[t][k]);
//            cv::waitKey();
        }//for(k)
    }//for(t)

    AQuA::opts.maxValueDat1 = mmax;
    AQuA::opts.minValueDat1 = mmin;
    normalizedParameter = static_cast<float>(mmax -mmin);

//    for (int t = 0; t < T; ++t) {
//        for (int k = 0; k < L; ++k) {
//            for (int i = 0; i < H; ++i) {
//                for (int j = 0; j < W; ++j) {
//                    frame[t][k].at<uchar>(i,j) -= static_cast<uchar>(AQuA::opts.minValueDat1);
//                }//for(i)
//            }//for(j)
//        }//for(k)
//    }//for(t)


    for (int t = 0; t < T; ++t) {
        for (int k = 0; k < L; ++k) {
            frame[t][k].convertTo(frame[t][k], CV_32F);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
//                    std::cout<< frame[t][k].at<float>(i,j);
                    frame[t][k].at<float>(i,j) -= static_cast<float>(AQuA::opts.minValueDat1);
//                    if (frame[t][k].at<float>(i,j)<0) {
//                        frame[t][k].at<float>(i, j) = 0;
//                    }
                    frame[t][k].at<float>(i,j) /= normalizedParameter;
                    std::cout<< frame[t][k].at<float>(i,j)<< "  ";
                }//for(i)
            }//for(j)
//            std::cout << "slice " << k + 1 << " of " << L << " frame " << t + 1 << " of " << T << std::endl;
//            cv::imshow("length", frame[t][k]);
//            cv::waitKey();
            std::cout<<std::endl;
        }//for(k)
        std::cout<<std::endl;
    }//for(t)

    AQuA::opts.sz[0] = H;
    AQuA::opts.sz[1] = W;
    AQuA::opts.sz[2] = L;
    AQuA::opts.sz[3] = T;
    AQuA::opts.BitDepth = BitDepth;

}//loadData()



int main(){
    AQuA::Init();
    loadData();
    return 0;
}

