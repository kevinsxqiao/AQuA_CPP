#include "data/data.h"
#include "preProcessRun/preProcessRun.h"

std::vector<cv::Mat> fit_F0_var_test(const std::vector<cv::Mat>& F0ProOrg, const std::vector<cv::Mat>& varMapOrg, int dist) {
    std::vector<cv::Mat> varMapOut(std::vector<cv::Mat>(L));
//        std::vector<cv::Mat> F0Pro (std::vector<cv::Mat>(L));
//        std::vector<cv::Mat> varMap (std::vector<cv::Mat>(L));
    std::vector<float> F0Pro;
    std::vector<float> varMap;
    std::vector<float> select1;
//        std::vector<cv::Mat> select1 (std::vector<cv::Mat>(L));
    std::vector<cv::Mat> select2(std::vector<cv::Mat>(L));

//#pragma omp parallel for
    for (int k = 0; k < L; ++k) {
//            F0Pro[k] = cv::Mat(H-2*dist,W-2*dist,CV_32F);
        varMapOut[k] = cv::Mat(H, W, CV_32F);
//            varMap[k] = cv::Mat(H-2*dist,W-2*dist,CV_32F);
    }

    //remove values close to 0
//#pragma omp parallel for collapse(3)
//        for (int k = 0; k < L; ++k) {
//            cv::Rect roi(dist,dist,H-dist,W-dist);
//            F0Pro[k] = F0ProOrg[k](roi).clone();
//            varMap[k] = varMapOrg[k](roi).clone();
//            for (int i = dist; i < H-dist; ++i) {
//                for (int j = dist; j < W-dist; ++j) {
//                    if (varMapOrg[k].at<float>(i,j) < 1e-8) {
//                        F0Pro[k].at<float>(i, j) = NAN;
//                        varMap[k].at<float>(i, j) = NAN;
//                    }
//                }
//            }
//        }
    for (int k = 0; k < L; ++k) {
        for (int i = dist; i < H - dist; ++i) {
            for (int j = dist; j < W - dist; ++j) {
                if (varMapOrg[k].at<float>(i, j) > 1e-8) {
                    F0Pro.push_back(F0ProOrg[k].at<float>(i, j));
                    varMap.push_back(varMapOrg[k].at<float>(i, j));
                }
            }
        }
    }
//        std::cout<<"F0Pro: "<<std::endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                std::cout<<F0Pro[0].at<float>(i,j)<<" ";
//            }
//            std::cout<<std::endl;
//        }
//        if isempty(F0Pro)
//           % if no valid variance
//        varMapOut = ones(size(F0ProOrg)) * 1e-8;
//        return;
//        end

    //downsample --- the number of original pixels is too large.

//        AQuA::writeDataToMatFile(F0Pro, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test_F0.mat");
//        cv::Point *maxLoc;
    float max_value = *std::max_element(F0Pro.begin(), F0Pro.end());
    float min_value = *std::min_element(F0Pro.begin(), F0Pro.end());

    if (dist == 1) {
        std::cout << "min max value:" << min_value << " " << max_value << std::endl;
//            std::cout<<"F0Pro: "<<std::endl;
//            for (int k = 0; k < 5; ++k) {
//                for (int i = 0; i < 20; ++i) {
//                    for (int j = 0; j < 20; ++j) {
//                        std::cout<<F0Pro[0].at<float>(i,j)<<" ";
//                    }
//                    std::cout<<std::endl;
//                }
//            }
    }
}


int main(){
    AQuA::Init();
    std::vector<std::vector<cv::Mat>> dataOrg = AQuA::loadData();
    auto start = std::chrono::high_resolution_clock::now();
    bool*** evtSpatialMask = AQuA::createEvtSpatialMask();
    std::vector<std::vector<cv::Mat>> dataNew = AQuA::baselineRemoveAndNoiseEstimation(dataOrg, evtSpatialMask);
    AQuA::writeDataToMatFile(dataNew, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test2.mat");
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    AQuA::release3dMatrix_bool(evtSpatialMask,H,W);
    return 0;
}