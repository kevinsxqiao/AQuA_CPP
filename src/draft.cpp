#include "data/data.h"
#include "actRun/actRun.h"
#include "preProcessRun/preProcessRun.h"
#include "data/data.h"
#include "phaseRun/phaseRun.h"
void func(){
    static vector<float> a;
    if (!a.empty()){
        cout<<1;
    }
    a.emplace_back(5);
}

int main(){
//    vector<vector<cv::Mat>> mat(500,vector<cv::Mat>(500));
//    for (int i = 0; i < 500; ++i) {
//        for (int j = 0; j < 500; ++j) {
//            mat[i][j] = cv::Mat::ones(5,5,CV_8U);
//            cv::multiply(mat[i][j],2,mat[i][j]);
//        }
//    }
//    vector<int> test ={2,3,4,5,6,7,8};
//    auto ite = find(test.begin(), test.end(),3);
//    int index = distance(test.begin(),ite);
//    cout<<test[index]<<" ";
    for (int i = 0; i < 5; ++i) {
        func();
    }
//    vector<int> curIt;
//    for (int ii_ite = 0; ii_ite < ihwOrg.size(); ++ii_ite) {
//        auto ite = find(ihwOrg.begin(), ihwOrg.end(),ihw[ii_ihw]);
//        if (ite != ihw.end()){
//            int index = distance(ihw.begin(),ite);
//            curIt.emplace_back(it[index]);
//        }
//    }//for(ii_ite)
//    vector<vector<cv::Mat>> select(200,vector<cv::Mat>(50));
//    for (int i = 0; i < 200; ++i) {
//        for (int j = 0; j < 50; ++j) {
//            select[i][j] = cv::Mat::zeros(95,95,CV_8U);
//        }
//    }
//    select[101][2].at<uchar>(89,6) = 1;
//    select[101][1].at<uchar>(89,5) = 1;
//
//    vector<vector<int>> res = AQuA::bw2Reg(select);
//    for (const auto& out:res) {
//        for(const auto& in:out){
//            cout<<in<<" ";
//        }
//        cout<<endl;
//    }
//
//    cout<<endl;
//    cout<<"correct ones:"<<endl;
//    cout<<AQuA::sub2ind(89,5,1,101,95,95,50)<< endl;
//    cout<<AQuA::sub2ind(89,6,2,101,95,95,50)<< endl;


//    cv::Mat mat = cv::Mat::zeros(5,5,CV_8U);
//    cv::Mat mask = cv::Mat::zeros(5,5,CV_8U);
//    mask.at<uchar>(0,1) = 1;
//    mask.at<uchar>(1,2) = 1;
//    mask.at<uchar>(2,3) = 1;
//    cv::Mat newmat = mat(mask);
//    for (int i = 0; i < newmat.total(); ++i) {
//        cout << static_cast<int>(newmat.at<uchar>(i)) << " ";
//    }
//    cout << endl;



//    mat[0][0].at<uchar>(0, 1) = 1;

//    mat[4][4].at<uchar>(4, 4) = 1;
//    for (int i = 0; i < mat[0][0].rows; ++i) {
//        for (int j = 0; j < mat[0][0].cols; ++j) {
//            std::cout<<static_cast<int>(mat[0][0].at<uchar>(i,j))<<"  ";
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<<endl;
//    cv::Mat mask1 = ((mat == 3) | (mat ==2));
//    vector<vector<AQuA::Point_struct>> res = AQuA::bw2Reg(mat);
//    mat.setTo(INFINITY,mask1);
//    cout<<res.size()<<" ";
//    cout<<res[0].size()<<" ";
//    cv::Mat dst = AQuA::myResize(mat,2,2);
//    for (int i = 0; i < dst.rows; ++i) {
//        for (int j = 0; j < dst.cols; ++j) {
//            std::cout<<dst.at<float>(i,j)<<"  ";
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<<mat.type()<<std::endl;
//    mat.at<uchar>(1,1) = 0;
//    mat.at<uchar>(2,2) = 0;
//    int flag = 0;
//    for (int i = 0; i < 5; ++i) {
//        for (int j = 0; j < 5; ++j) {
//            std::cout<<static_cast<int>(mat.at<uchar>(i,j))<<"  ";
//            if (mat.at<uchar>(i,j) !=1){
//                flag++;
//            }
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<<flag<<std::endl;
//    std::cout<<"size of 8U: "<<sizeof (uchar)<<std::endl;
//    std::cout<<"size of 8S: "<<sizeof (char)<<std::endl;
//    std::cout<<"size of 32S: "<<sizeof (int)<<std::endl;

    return 0;
}


