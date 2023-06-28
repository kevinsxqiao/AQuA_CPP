#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <mat.h>
#include <vector>
#include <boost/math/distributions/normal.hpp>

float obtainBias(){
    MATFile *pmatFile;
    mxArray *pMxbiasMatrix; //215*55 double
    mxArray *pMxcuts; //1*55 double
    mxArray *pMxwindowSizes; //1*215 double
    double *biasMatrix;
    double *cuts;
    double *windowSizes;

    std::cout<< "--------loading cfg F0_biasMatrix--------"<<std::endl;
    const char *filename = "../cfg/F0_biasMatrix.mat";
//    const char * filename1 = "C:/Users/Kevin Qiao/Desktop/AQuA_CPP/cfg/F0_biasMatrix.mat";
    pmatFile = matOpen(filename, "r");
    if (pmatFile == nullptr) {
        std::cout<< "--------error opening cfg--------"<<std::endl;
        std::exit(-1);
    }

    pMxbiasMatrix = matGetVariable(pmatFile, "biasMatrix");
    if (pMxbiasMatrix == nullptr) {
        std::cout<< "--------error reading variable \"biasMatrix\" from file--------"<<std::endl;
        std::exit(-1);
    }
    pMxcuts = matGetVariable(pmatFile, "cuts");
    if (pMxcuts == nullptr) {
        std::cout<< "--------error reading variable \"cuts\" from file--------"<<std::endl;
        std::exit(-1);
    }
    pMxwindowSizes = matGetVariable(pmatFile, "windowSizes");
    if (pMxwindowSizes == nullptr) {
        std::cout<< "--------error reading variable \"windowSizes\" from file--------"<<std::endl;
        std::exit(-1);
    }

    biasMatrix = mxGetPr(pMxbiasMatrix);
    if (biasMatrix == nullptr) {
        std::cout<< "--------error reading data from variable \"biasMatrix\"-------"<<std::endl;
        std::exit(-1);
    }
    cuts = mxGetPr(pMxcuts);
    if (cuts == nullptr) {
        std::cout<< "--------error reading data from variable \"cuts\"-------"<<std::endl;
        std::exit(-1);
    }
    windowSizes = mxGetPr(pMxwindowSizes);
    if (windowSizes == nullptr) {
        std::cout<< "--------error reading data from variable \"windowSizes\"-------"<<std::endl;
        std::exit(-1);
    }

    const mwSize *dims_biasMatrix = mxGetDimensions(pMxbiasMatrix);
    const mwSize *dims_cuts = mxGetDimensions(pMxcuts);
    const mwSize *dims_windowSizes = mxGetDimensions(pMxwindowSizes);

    float bias=0;

//        H = dims[0];
//        std::cout<<"original size: "<< std::endl;
//        frame[t][k].at<float>(i,j) = static_cast<float>(pdata[j_src*H + i_src + k*H*W + t*H*W*L]);
    if (25 > 200){
        bias = 0;
        return bias;
    }

    //linear interpolate value if cannot find it
    int idx0=-1, idy0=-1, idx1=-1, idy1=-1;
    for (int i = dims_windowSizes[1] - 1; i >= 0; --i) {
        if (windowSizes[i] <= 25){
            idx0 = i;
            break;
        }
    }
    for (int i = dims_cuts[1] - 1; i >= 0; --i) {
        if (cuts[i] <= 200){
            idy0 = i;
            break;
        }
    }
    if ((idx0==dims_windowSizes[1]-1) || (windowSizes[idx0]==25)){
        idx1 = idx0;
    } else{
        idx1 = idx0 + 1;
    }
    if ((idy0==dims_cuts[1]-1) || (cuts[idy0]==200)){
        idy1 = idy0;
    } else{
        idy1 = idy0 + 1;
    }

    float bias0=0, bias1=0;
    for (int i = idx0; i <= idx1 ; ++i) {
        if (!std::isnan(biasMatrix[i + idy0 * dims_biasMatrix[0]])){
            bias0 += static_cast<float>(biasMatrix[i + idy0 * dims_biasMatrix[0]]);
        }
        if (!std::isnan(biasMatrix[i + idy1 * dims_biasMatrix[0]])){
            bias1 += static_cast<float>(biasMatrix[i + idy1 * dims_biasMatrix[0]]);
        }
    }
    bias0 /= static_cast<float>(idx1 - idx0 + 1);
    bias1 /= static_cast<float>(idx1 - idx0 + 1);
    int cut0 = static_cast<int>(cuts[idy0]);
    int cut1 = static_cast<int>(cuts[idy0+1]);

    if (std::isnan(bias0)){
        bias = bias1;
    } else{
        bias = bias0 + (bias1 - bias0) / (cut1 - cut0) * (200 - cut0);
    }

    return bias;
}//obtainBias()

int main(){
    float bias;
    bias = obtainBias();
    return 0;
}

//int main() {
//    cv::Mat mat = (cv::Mat_<double>(3, 3) << 1, 2, 3,
//            4, 5, 6,
//            7, 8, 9);
//    std::vector<cv::Mat> matrix(2);
//    for (int i = 0; i < 2; ++i) {
//        matrix[i] = mat.clone();
//    }
//    std::vector<std::vector<cv::Mat>> matri(2,std::vector<cv::Mat> (2));
//    for (int i = 0; i < 2; ++i) {
//        for (int j = 0; j < 2; ++j) {
//            matri[i][j] = matrix[j].clone();
//        }
//    }
////    std::cout<< matri[0][1].at<double>(2,2);
//    std::cout<< matri.size();
//    std::cout<< matri[1].size();
//    std::cout<< matri[1][1].size();
////    double minVal, maxVal;
////    cv::Point minLoc, maxLoc;
////
////    cv::minMaxLoc(mat, &minVal, &maxVal, &minLoc, &maxLoc);
////
////    std::cout << "Min value: " << minVal << " at (" << minLoc.x << ", " << minLoc.y << ")\n";
////    std::cout << "Max value: " << maxVal << " at (" << maxLoc.x << ", " << maxLoc.y << ")\n";
//
//    return 0;
//}

