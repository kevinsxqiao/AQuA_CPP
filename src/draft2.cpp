vector<cv::Mat> load3D(const char* fileName, const char* varName) {
    MATFile *pmatFile;
    mxArray *pMxArray;

    cout<< "--------loading data--------"<<endl;
    pmatFile = matOpen(fileName, "r");
    if (pmatFile == nullptr) {
        cout<< "--------error opening file--------"<<endl;
        exit(-1);
    }

    pMxArray = matGetVariable(pmatFile, varName);
    if (pMxArray == nullptr) {
        cout<< "--------error reading variable from file--------"<<endl;
        exit(-1);
    }

    void* pdata = mxGetData(pMxArray);
    mxClassID classID = mxGetClassID(pMxArray);

    if (pdata == nullptr) {
        cout<< "--------error reading data from variable-------"<<endl;
        exit(-1);
    }

    const mwSize *dims = mxGetDimensions(pMxArray);

    vector<cv::Mat> frame(dims[2]);

    if (classID == mxSINGLE_CLASS){
        for (int k = 0; k < dims[2]; ++k) {
            frame[k] = cv::Mat(dims[0],dims[1],CV_32F);
            for (int i = 0; i < dims[0]; ++i) {
                for (int j = 0; j < dims[1]; ++j){
                    frame[k].at<float>(i,j) = static_cast<float*>(pdata)[sub2ind(i,j,k,dims[0],dims[1])];
                }//for(j)
            }//for(i)
        }//for(k)
    } else if(classID == mxDOUBLE_CLASS){
        for (int k = 0; k < dims[2]; ++k) {
            frame[k] = cv::Mat(dims[0],dims[1],CV_32F);
            for (int i = 0; i < dims[0]; ++i) {
                for (int j = 0; j < dims[1]; ++j){
                    frame[k].at<float>(i,j) = static_cast<double*>(pdata)[sub2ind(i,j,k,dims[0],dims[1])];
                }//for(j)
            }//for(i)
        }//for(k)
    } else{
        cout << "Unhandled data type." << endl;
    }

    cout<<"height of image:"<< dims[0] << endl;
    cout<<"width of image:"<< dims[1] << endl;
    cout<<"length of image:"<< dims[2] << endl;
    cout<<"--------data loaded--------"<<endl;

    //release MAT pointer
    if (pMxArray != nullptr) {
        mxDestroyArray(pMxArray);
    }

    if (pmatFile != nullptr) {
        matClose(pmatFile);
    }

    return frame;
}//load4D()