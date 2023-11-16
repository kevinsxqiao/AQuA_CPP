#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <mat.h>
#include <vector>
#include <boost/math/distributions/normal.hpp>
#include "data/data.h"
#include <omp.h>
#include "phaseRun/phaseRun.h"
#include "actRun/actRun.h"
#include <chrono>





int main() {
    vector<vector<cv::Mat>> sel_mat = AQuA::load4D_8U("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/sel_mat.mat", "selectMap");
//    vector<vector<cv::Mat>> sel_cpp = AQuA::load4D_8U("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/sel_cpp.mat", "myVar");
//    vector<vector<cv::Mat>> dF1 = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/phaseRun.mat", "dF1");
//    AQuA::writeDataToMatFile(dF1, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test/dF_test.mat");
//    vector<vector<cv::Mat>> dF1_test = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/dF_test.mat", "myVar");
//    vector<vector<cv::Mat>> dF_cpp = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/dFRe0.mat", "myVar");
//    vector<vector<cv::Mat>> val_cpp = AQuA::load4D_8U("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/val0.mat", "myVar");
//    vector<vector<cv::Mat>> dF_mat = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/dFRe0_mat.mat", "mat11");
//    vector<vector<cv::Mat>> val_mat = AQuA::load4D_8U("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/val0_mat.mat", "mat1");
//            float a = 5.0;
//            uchar b = static_cast<uchar>(ceil((a/2.0)));
//            cout<<static_cast<float>(b);
//        cv::Mat test2D = cv::Mat::ones(2,3,CV_32F);
//        test2D.at<float>(0,1) = 0;
//        AQuA::writeDataToMatFile(test2D, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test/test2D.mat");
//    for (int t = 0; t < dF_cpp.size(); ++t) {
//        for (int k = 0; k < dF_cpp[0].size(); ++k) {
//            cv::Mat dst_dF;
//            cv::Mat dst_val;
////            cv::bitwise_xor(dF_cpp[t][k],dF_mat[t][k],dst_dF);
//            cv::bitwise_xor(val_cpp[t][k],val_mat[t][k],dst_val);
////            if (cv::countNonZero(dst_dF) >0){
////                cout<<"dF: "<<"t: "<<t<<" k:"<<k<<endl;
////            }
//            if (cv::countNonZero(dst_val) >0){
//                cout<<"val: "<<"t: "<<t<<" k:"<<k<<endl;
//            }
//        }
//    }
    float scaleRatio =8;
    vector<vector<int>> curRegions = AQuA::bw2Reg(sel_mat);
    curRegions.erase(remove_if(curRegions.begin(), curRegions.end(),
                               [&](const auto& region) {
                                   return region.size() <= (20 / pow(8, 2) * 5 / 3);
                               }), curRegions.end()); //move all the wanted elements to the front, and iterator points at next pos of the element wanted, then delete


    for (int ii_cur = 0; ii_cur < curRegions.size(); ++ii_cur) { //--i
        unordered_set<int> ihw_temp;
        vector<int> pix;
        vector<int> ih;
        vector<int> iw;
        vector<int> il;
        vector<int> it;
        for (int ii_ite = 0; ii_ite < curRegions[ii_cur].size(); ++ii_ite) {
            pix.emplace_back(curRegions[ii_cur][ii_ite]);
        }
        sort(pix.begin(), pix.end());

        for (int ii_ite = 0; ii_ite < pix.size(); ++ii_ite) {
            AQuA::Point_struct pix_temp = AQuA::ind2sub(pix[ii_ite],24,24,50);
            ih.emplace_back(pix_temp.i);
            iw.emplace_back(pix_temp.j);
            il.emplace_back(pix_temp.k);
            it.emplace_back(pix_temp.t);
            ihw_temp.insert(AQuA::sub2ind(ih[ii_ite],iw[ii_ite],il[ii_ite],24,24));
        }
        vector<int> ihw(ihw_temp.begin(), ihw_temp.end());
        sort(ihw.begin(), ihw.end());
        int dur = *max_element(it.begin(),it.end()) - *min_element(it.begin(), it.end()) + 1;
        int arLabel_hs = ih[0] * 8;
        int arLabel_he = min(static_cast<float>(190),(ih[0]+1)*scaleRatio - 1);
        int arLabel_ws = iw[0] * 8;
        int arLabel_we = min(static_cast<float>(190),(iw[0]+1)*scaleRatio - 1);
        vector<int> arLabel;
//        for (int i = arLabel_hs; i <= arLabel_he; ++i) {
//            for (int j = arLabel_ws; j <= arLabel_we; ++j) {
//                if (activeMap[it[0]][il[0]].at<ushort>(i, j) != 0) {
//                    arLabel.emplace_back(static_cast<int>(activeMap[it[0]][il[0]].at<ushort>(i, j)));
//                }
//            }
//        }
        int arLabel_val = *min_element(arLabel.begin(), arLabel.end());
//    vector<vector<int>> cur_cpp = AQuA::bw2Reg(sel_cpp);
    return 0;
}
}
