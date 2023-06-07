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
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<cv::Mat>> data = AQuA::loadData();
    AQuA::baselineRemoveAndNoiseEstimation(data);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "used time: " << duration/1000 << " seconds" << std::endl;
    return 1;
}