int main(){
    AQuA::Init();
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> data1 = AQuA::loadData();
    AQuA::regCrossCorrelation(data1);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_mat = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "time: " << duration_mat << " milliseconds" << std::endl;
//    std::cout<< AQuA::opts.minValueDat1 << "  "<< AQuA::opts.maxValueDat1;
//    std::cout<< data1[0][0].at<float>(0,0);
    return 0;
}


int main(){
    AQuA::Init();
    H=2;
    W=2;
    L=2;
    float a[3][3][3] = {
            {{1,2,3}, {4,5,0}, {0,1,2}},
            {{3,4,5}, {0,1,2}, {2,3,4}},
            {{5,0,1}, {2,3,4}, {4,5,0}}
    };

    float b[3][3][3] = {
            {{5, 4, 3}, {2, 1, 0}, {5, 4, 3}},
            {{3, 2, 1}, {0, 5, 4}, {1, 0, 5}},
            {{1, 0, 5}, {4, 3, 2}, {1, 0, 5}}
    };
    float*** a_add = AQuA::create3dMatrix(3,3,3);
    float*** b_add = AQuA::create3dMatrix(3,3,3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                a_add[i][j][k]= a[i][j][k];
                b_add[i][j][k]= b[i][j][k];
//                std::cout<< a_add[i][j][k]<< " "<< b_add[i][j][k]<< " " ;
            }
        }
    }
    AQuA::dft(a_add, b_add);
    AQuA::release3dMatrix(a_add,3,3);
    AQuA::release3dMatrix(b_add,3,3);
    return 1;
}

int main(){
    AQuA::Init();
    std::vector<std::vector<cv::Mat>> dataOrg = AQuA::loadData();
    auto start = std::chrono::high_resolution_clock::now();
    bool*** evtSpatialMask = AQuA::createEvtSpatialMask();
    std::vector<std::vector<cv::Mat>> dataNew = AQuA::preProcessRun(dataOrg, evtSpatialMask);
    AQuA::writeDataToMatFile(dataNew, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test2.mat");
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    AQuA::release3dMatrix_bool(evtSpatialMask,H,W);
    return 0;
}

int main(){
    std::vector<std::vector<cv::Mat>> mat(5,std::vector<cv::Mat>(5));
    for (int t = 0; t < 5; ++t) {
        for (int k = 0; k < 5; ++k) {
            mat[t][k] = cv::Mat::zeros(5,5,CV_32S);
        }
    }

    mat[0][0].at<int>(0,0) = 1;
//    mat[0][0].at<int>(1,0) = 1;
//    mat[0][1].at<int>(1,2) = 1;

    mat[0][1].at<int>(4,4) = 1;
    mat[0][1].at<int>(3,4) = 1;
//    mat[0][2].at<int>(1,2) = 1;
//    mat[0][2].at<int>(1,0) = 1;

    mat[4][0].at<int>(1,1) = 1;
    std::vector<std::vector<AQuA::Point_struct>> poi = AQuA::bw2Reg(mat);
    std::cout<< poi.size()<<std::endl;
    std::cout<< poi[0].size()<<std::endl;
    std::cout<< poi[1].size()<<std::endl;
    std::cout<< poi[1][0].t<< poi[1][0].k<< poi[1][0].i<< poi[1][0].j<<std::endl;


    return 0;
}