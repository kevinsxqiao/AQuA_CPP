template<typename T>
std::vector<std::vector<cv::Mat>> loadData(T* pdata, const mwSize* dims, mwSize numDims) {
    std::vector<std::vector<cv::Mat>> frame;

    if (numDims == 4){
        H = dims[0];
        W = dims[1];
        L = dims[2];
        T = dims[3];
    }
    if (numDims == 3){
        H = dims[0];
        W = dims[1];
        L = dims[2];
    }

    frame.resize(T,std::vector<cv::Mat>(L));
    for (int t = 0; t < T; ++t) {
        for (int k = 0; k < L; ++k) {
            frame[t][k] = cv::Mat(H,W,CV_32F);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j){
                    frame[t][k].at<T>(i,j) = pdata[sub2ind(i,j,k,t,H,W,L)];
                }
            }
        }
    }

    return frame;
}

std::vector<std::vector<cv::Mat>> load4DData_clean(const char* fileName, const char* varName) {
    // 读取MAT文件，检查数据类型等代码...

    mxClassID classID = mxGetClassID(pMxArray);
    void* pdata = mxGetData(pMxArray);

    std::vector<std::vector<cv::Mat>> frame;
    if (classID == mxSINGLE_CLASS) {
        frame = loadData(static_cast<float*>(pdata), dims, numDims);
    } else if (classID == mxDOUBLE_CLASS) {
        frame = loadData(static_cast<double*>(pdata), dims, numDims);
    } else {
        std::cout << "Unhandled data type." << std::endl;
    }

    // 释放MAT指针等代码...

    return frame;
}
