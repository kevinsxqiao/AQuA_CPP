//
// Created by Kevin Qiao on 3/28/2023.
//

#include "data.h"

class path;

class path;

using namespace std;

namespace AQuA{

    /*
     * define data structure
     */
//    struct rawDataSize_struct rawDataSize;
    struct preSetting_struct preSetting;
    struct opts_struct opts;



    //    judge if registration and bleach have been executed; ---- true = no ; false = both executed
    bool isDefault() {
        return (preSetting.registrateCorrect == preSetting.registrateCorrect_default
                && preSetting.bleachCorrect == preSetting.bleachCorrect_default);
    }// isDefault()


    void preSettingInit(){
        preSetting.registrateCorrect = preSetting.registrateCorrect_default;
        preSetting.bleachCorrect = preSetting.bleachCorrect_default;
        cout<< "--------preSetting initialized--------"<<endl;
    }// preSettingInit()


    void optsInit(){
        opts.fileName1 = "C:/Users/Kevin Qiao/Desktop/AQuA_data/Test_global_local_3D.mat";
        opts.alreadyPreprocess = false;
        opts.alreadyBleachCorrect = false;
        opts.movAvgWin = 25;
        opts.cut = 200;
        opts.smoXY = 0.5;
        opts.BitDepth = -1;
        opts.regMaskGap = 5;
        opts.singleChannel = true;
        opts.registrateCorrect = 0;
        opts.bleachCorrect = 0;
        opts.medSmo = 1;
        opts.cut = 200;
        opts.maxdF1 = 1;
        opts.thrARScl = 3;
        opts.minSize = 20;
        opts.maxSize = INFINITY;
        opts.minDur = 5;
        opts.circularityThr = 0;
        opts.spaMergeDist = 0;
        opts.compress = 0;
        opts.needTemp = true;
        opts.step = 0.5;
        opts.sigThr = 3.5;
        opts.maxDelay = 0.6;
        opts.seedSzRatio = 3.5;
        opts.needRefine = false;
        opts.needGrow = false;
        opts.maxSpaScale = 7;
        opts.minSpaScale = 3;
        opts.TPatch = 20;

        cout<< "--------opts initialized--------"<<endl;
    }// optsInit()


//    void rawDataSizeInit(){
//        rawDataSize={0};
//        cout<< "--------rawDataSize initialized--------"<<endl;
//    }


    void Init(){
        preSettingInit();
        optsInit();
//        rawDataSizeInit();
        cout<<endl;
    }


    float medianFunc(float* array, int size){
        float median;
        sort(array, array+size);
        if(size % 2 == 0){
            median = ( array[size/2] + array[size/2 - 1] ) / 2;
        }
        else{
            median = array[size/ 2];
        }
//        cout<<"median value of reference frame: "<< median<<endl;
        return median;
    }//medianFunc()


    vector<vector<cv::Mat>> loadData() {
        MATFile *pmatFile;
        mxArray *pMxArray;
        double *pdata;
        int bdCrop = AQuA::opts.regMaskGap;
        double min,max=0;
        double mmin=255, mmax=0;
        int BitDepth = -1;
        float normalizedParameter;

        cout<< "--------loading data--------"<<endl;
//        cout<<'\r'<< 0 << "% "<<flush;
//        const char *filename = "C:/Users/Kevin Qiao/Desktop/AQuA_data/Test_global_local_3D.mat";
        pmatFile = matOpen(opts.fileName1, "r");
        if (pmatFile == nullptr) {
            cout<< "--------error opening file--------"<<endl;
            exit(-1);
        }

        pMxArray = matGetVariable(pmatFile, "synthetic");
        if (pMxArray == nullptr) {
            cout<< "--------error reading variable from file--------"<<endl;
            exit(-1);
        }

        pdata = mxGetPr(pMxArray);
        if (pdata == nullptr) {
            cout<< "--------error reading data from variable-------"<<endl;
            exit(-1);
        }

        const mwSize *dims = mxGetDimensions(pMxArray);
        int H = dims[0];
        int W = dims[1];
        int L = dims[2];
        int T = dims[3];
//        T = 10;
//        cout<<"original size: "<< endl;
//        cout<<"height of image:"<< H << endl;
//        cout<<"width of image:"<< W << endl;
//        cout<<"length of image:"<< L << endl;
//        cout<<"time frames of image:"<< T << endl;

        vector<vector<cv::Mat>> frame(T);
        for (int t = 0; t < T; ++t) {
//            cout<<'\r'<< 100*t/T/2 <<"% "<<flush;
            for (int k = 0; k < L; ++k) {
                frame[t].emplace_back(H- 2*bdCrop,W- 2*bdCrop,CV_32F);
//#pragma omp parallel for collapse(2)
                for (int i = 0, i_src = bdCrop; i < H - 2*bdCrop; ++i, ++ i_src) {
                    for (int j = 0, j_src= bdCrop; j < W - 2*bdCrop; ++j, ++j_src) {
                        frame[t][k].at<float>(i,j) = static_cast<float>(pdata[j_src*H + i_src + k*H*W + t*H*W*L]);
                    }//for(j)
                }//for(i)
                cv::minMaxLoc(frame[t][k], &min, &max);
                if (min < mmin){
                    mmin = min;
                }
                if(max > mmax){
                    mmax = max;
                }
            }//for(k)
        }//for(t)

        H = dims[0] - 2*bdCrop;
        W = dims[1] - 2*bdCrop;


        //release MAT pointer
        if (pMxArray != nullptr) {
            mxDestroyArray(pMxArray);
        }

        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }

        AQuA::opts.maxValueDat1 = mmax;
        AQuA::opts.minValueDat1 = mmin;
//        cout<<"minValue: "<< opts.minValueDat1 << "  maxValue: "<<opts.maxValueDat1<<endl;
        normalizedParameter = static_cast<float>(mmax -mmin);

//#pragma omp parallel for collapse(4)
        for (int t = 0; t < T; ++t) {
//            cout<<'\r'<< 50+ (100*t/T/2) <<"%"<<flush;
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        frame[t][k].at<float>(i,j) = (frame[t][k].at<float>(i,j) - static_cast<float>(AQuA::opts.minValueDat1)) / normalizedParameter;
//                    cout<< frame[t][k].at<float>(i,j)<< "  ";// display pixel value
                    }//for(i)
                }//for(j)
            }//for(k)
        }//for(t)
//        cout << "\r" << 100 << "%" << flush<< endl;
        AQuA::opts.sz[0] = H;
        AQuA::opts.sz[1] = W;
        AQuA::opts.sz[2] = L;
        AQuA::opts.sz[3] = T;
        AQuA::opts.BitDepth = BitDepth;

//        cout<<"data: "<<endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<< frame[0][0].at<float>(i,j)<< "  ";
//            }
//            cout<<endl;
//        }
    cout<<"after cropping: "<< endl;
    cout<<"height of image:"<< H << endl;
    cout<<"width of image:"<< W << endl;
    cout<<"length of image:"<< L << endl;
    cout<<"time frames of image:"<< T << endl;
    cout<<"--------data loaded--------"<<endl;
    opts.data1_org = frame;

//    for (int t = 0; t < T; ++t) {
//        vector<cv::Mat> dat(L);
//        for (int k = 0; k < L; ++k) {
//            dat[t].emplace_back(frame[t][k].clone());
//        }
//        opts.data1_org.emplace_back(dat);
//    }

    return frame;
    }//loadData()


    vector<vector<cv::Mat>> load4D(const char* fileName, const char* varName) {
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


        vector<vector<cv::Mat>> frame(dims[3],vector<cv::Mat>(dims[2]));

        if (classID == mxSINGLE_CLASS){
            for (int t = 0; t < dims[3]; ++t) {
                for (int k = 0; k < dims[2]; ++k) {
                    frame[t][k] = cv::Mat(dims[0],dims[1],CV_32F);
                    for (int i = 0; i < dims[0]; ++i) {
                        for (int j = 0; j < dims[1]; ++j){
                            frame[t][k].at<float>(i,j) = static_cast<float*>(pdata)[sub2ind(i,j,k,t,dims[0],dims[1],dims[2])];
                        }//for(j)
                    }//for(i)
                }//for(k)
            }//for(t)
        } else if(classID == mxDOUBLE_CLASS){
            for (int t = 0; t < dims[3]; ++t) {
                for (int k = 0; k < dims[2]; ++k) {
                    frame[t][k] = cv::Mat(dims[0],dims[1],CV_32F);
                    for (int i = 0; i < dims[0]; ++i) {
                        for (int j = 0; j < dims[1]; ++j){
                            frame[t][k].at<float>(i,j) = static_cast<double*>(pdata)[sub2ind(i,j,k,t,dims[0],dims[1],dims[2])];
                        }//for(j)
                    }//for(i)
                }//for(k)
            }//for(t)
        } else{
            cout << "Unhandled data type." << endl;
        }

//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<frame[0][0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }

        cout<<"height of image:"<< dims[0] << endl;
        cout<<"width of image:"<< dims[1] << endl;
        cout<<"length of image:"<< dims[2] << endl;
        cout<<"time frames of image:"<< dims[3] << endl;
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
    }//load3D()
    
    
    cv::Mat load2D(const char* fileName, const char* varName) {
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

        cv::Mat frame = cv::Mat(dims[0],dims[1],CV_32F);

        if (classID == mxSINGLE_CLASS){
            for (int i = 0; i < dims[0]; ++i) {
                for (int j = 0; j < dims[1]; ++j){
                    frame.at<float>(i,j) = static_cast<float*>(pdata)[i+j*dims[0]];
                }//for(j)
            }//for(i)
        } else if(classID == mxDOUBLE_CLASS){
            for (int i = 0; i < dims[0]; ++i) {
                for (int j = 0; j < dims[1]; ++j){
                    frame.at<float>(i,j) = static_cast<double*>(pdata)[i+j*dims[0]];
                }//for(j)
            }//for(i)
        } else{
            cout << "Unhandled data type." << endl;
        }

        cout<<"height of image:"<< dims[0] << endl;
        cout<<"width of image:"<< dims[1] << endl;

        //release MAT pointer
        if (pMxArray != nullptr) {
            mxDestroyArray(pMxArray);
        }

        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }


        cout<<"--------data loaded--------"<<endl;
        return frame;
    }//load2D()


    vector<vector<int>> loadCell_int(const char* fileName, const char* varName){
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

        mwSize numCells = mxGetNumberOfElements(pMxArray);
        mxArray* pCellElement = mxGetCell(pMxArray, 0);
        mxClassID classID = mxGetClassID(pCellElement);//get data type in cell

        vector<vector<int>> frame(numCells);

        if (classID == mxSINGLE_CLASS){
            for (mwIndex i = 0; i < numCells; i++) {
                pCellElement = mxGetCell(pMxArray, i);
                if (pCellElement != NULL) {
                    mwSize numElements = mxGetNumberOfElements(pCellElement);
                    float* pdata = static_cast<float*>(mxGetData(pCellElement));

                    frame[i].resize(numElements);
                    transform(pdata, pdata + numElements, frame[i].begin(), [](float val) { return static_cast<int>(val); });
                }
            }
        } else if(classID == mxDOUBLE_CLASS){
            for (mwIndex i = 0; i < numCells; i++) {
                pCellElement = mxGetCell(pMxArray, i);
                if (pCellElement != NULL) {
                    mwSize numElements = mxGetNumberOfElements(pCellElement);
                    double* pdata = static_cast<double*>(mxGetData(pCellElement));

                    frame[i].resize(numElements);
                    transform(pdata, pdata + numElements, frame[i].begin(), [](double val) { return static_cast<int>(val); });
                }
            }
        } else{
            cout << "Unhandled data type." << endl;
        }
        cout<<"number of cells:"<< numCells << endl;
        //release MAT pointer
        if (pMxArray != nullptr) {
            mxDestroyArray(pMxArray);
        }

        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }
        cout<<"data loaded"<< endl;
        return frame;
    }//loadCell_int()


    vector<vector<double>> loadCell_double(const char* fileName, const char* varName){
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

        mwSize numCells = mxGetNumberOfElements(pMxArray);
        mxArray* pCellElement = mxGetCell(pMxArray, 0);
//        mxClassID classID = mxGetClassID(pCellElement);//get data type in cell

        vector<vector<double>> frame(numCells);
        for (mwIndex i = 0; i < numCells; i++) {
            pCellElement = mxGetCell(pMxArray, i);
            if (pCellElement != NULL) {
                mwSize numElements = mxGetNumberOfElements(pCellElement);
                double* pdata = static_cast<double*>(mxGetData(pCellElement));

                frame[i].resize(numElements);
                transform(pdata, pdata + numElements, frame[i].begin(), [](double val) { return (val); });
            }
        }
        cout<<"number of cells:"<< numCells << endl;
        //release MAT pointer
        if (pMxArray != nullptr) {
            mxDestroyArray(pMxArray);
        }

        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }
        cout<<"data loaded"<< endl;
        return frame;
    }//loadCell_double()


    vector<cv::Mat> loadCell_matrix(const char* fileName, const char* varName){
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

        mwSize numCells = mxGetNumberOfElements(pMxArray);
        mxArray* pCellElement = nullptr;

        vector<cv::Mat> matrix(numCells);

        for (mwIndex k = 0; k < numCells; k++) {
            pCellElement = mxGetCell(pMxArray, k);
            if (pCellElement != nullptr) {
                void *pdata = (mxGetData(pCellElement));
//                mxClassID classID = mxGetClassID(pCellElement);
                const mwSize *dims = mxGetDimensions(pCellElement);

                matrix[k] = cv::Mat(dims[0], dims[1], CV_32F);
//                std::memcpy(matrix[k].data, pdata, dims[0] * dims[1] * sizeof(float));
                for (int i = 0; i < dims[0]; ++i) {
                    for (int j = 0; j < dims[1]; ++j){
                        matrix[k].at<float>(i,j) = static_cast<double*>(pdata)[i+j*dims[0]];
                    }//for(j)
                }//for(i)
            }//if
        }//for(k)
        cout<<"number of cells:"<< numCells << endl;
        //release MAT pointer
        if (pMxArray != nullptr) {
            mxDestroyArray(pMxArray);
        }

        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }
        cout<<"data loaded"<< endl;
        return matrix;
    }//loadCell_matrix()


    mxArray* cvDataToMxArray(const cv::Mat& data) {
        // Calculate the size of the 2D matrix
        mwSize dims[2] = {static_cast<mwSize>(data.rows), static_cast<mwSize>(data.cols)};

        // Create a 2D mxArray
        mxArray* pMxArray = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);

        // Copy data from your vector to the mxArray
        float* ptr = reinterpret_cast<float*>(mxGetData(pMxArray));
        memcpy(ptr, data.data, data.rows * data.cols * sizeof(float));

        return pMxArray;
    }


    mxArray* cvDataToMxArray(const vector<vector<cv::Mat>>& data) {
        // Calculate the size of the 4D matrix
        mwSize dims[4] = {static_cast<mwSize>(data[0][0].cols), static_cast<mwSize>(data[0][0].rows),
                          static_cast<mwSize>(data[0].size()), static_cast<mwSize>(data.size())};

        // Create a 4D mxArray
        mxArray* pMxArray = mxCreateNumericArray(4, dims, mxUINT8_CLASS, mxREAL);

        // Copy data from your vector to the mxArray
        uchar* ptr = reinterpret_cast<uchar*>(mxGetUint8s(pMxArray));
        for (int t = 0; t < data.size(); ++t) {
            for (int k = 0; k < data[t].size(); ++k) {
                const cv::Mat& mat = data[t][k];
                for (int i = 0; i < mat.cols; ++i) {
                    for (int j = 0; j < mat.rows; ++j) {
                        *ptr = mat.at<uchar>(j, i);
                        ++ptr;
                    }
                }
            }
        }

        return pMxArray;
    }



    mxArray* cvDataToMxArray(const vector<cv::Mat>& data) {
        // Calculate the size of the 3D matrix
        mwSize dims[3] = {static_cast<mwSize>(data[0].rows), static_cast<mwSize>(data[0].cols), static_cast<mwSize>(data.size())};

        // Create a 3D mxArray
        mxArray* pMxArray = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);

        // Copy data from your vector to the mxArray
        float* ptr = reinterpret_cast<float*>(mxGetData(pMxArray));
            for (int k = 0; k < data.size(); ++k) {
                const cv::Mat& mat = data[k];
                memcpy(ptr, mat.data, mat.rows * mat.cols * sizeof(float));
                ptr += mat.rows * mat.cols;
            }
        return pMxArray;
    }


    void writeDataToMatFile(cv::Mat& data, const string& filename) {
        cout<<"--------start writing--------"<<endl;
        MATFile *pmatFile;

        // Open the mat file
        pmatFile = matOpen(filename.c_str(), "w");
        if (pmatFile == nullptr) {
            cout << "--------error opening file--------" << endl;
            exit(-1);
        }

        // Convert your data to a mxArray
        mxArray* pMxArray = cvDataToMxArray(data);

        // Write the variable to the mat file
        if (matPutVariable(pmatFile, "myVar", pMxArray) != 0) {
            cout << "--------error writing variable to file--------" << endl;
            exit(-1);
        }

        // Free the mxArray
        mxDestroyArray(pMxArray);

        // Close the mat file
        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }
        cout<<"--------finish writing--------"<<endl;
    }


    void writeDataToMatFile(vector<vector<cv::Mat>>& data, const string& filename) {
        cout<<"--------start writing--------"<<endl;
        MATFile *pmatFile;

        // Open the mat file
        pmatFile = matOpen(filename.c_str(), "w");
        if (pmatFile == nullptr) {
            cout << "--------error opening file--------" << endl;
            exit(-1);
        }

        // Convert your data to a mxArray
        mxArray* pMxArray = cvDataToMxArray(data);

        // Write the variable to the mat file
        if (matPutVariable(pmatFile, "myVar", pMxArray) != 0) {
            cout << "--------error writing variable to file--------" << endl;
            exit(-1);
        }

        // Free the mxArray
        mxDestroyArray(pMxArray);

        // Close the mat file
        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }
        cout<<"--------finish writing--------"<<endl;
    }


    void writeDataToMatFile(vector<cv::Mat>& data, const string& filename) {
        cout<<"--------start writing--------"<<endl;
        MATFile *pmatFile;

        // Open the mat file
        pmatFile = matOpen(filename.c_str(), "w");
        if (pmatFile == nullptr) {
            cout << "--------error opening file--------" << endl;
            exit(-1);
        }

        // Convert your data to a mxArray
        mxArray* pMxArray = cvDataToMxArray(data);

        // Write the variable to the mat file
        if (matPutVariable(pmatFile, "myVar", pMxArray) != 0) {
            cout << "--------error writing variable to file--------" << endl;
            exit(-1);
        }

        // Free the mxArray
        mxDestroyArray(pMxArray);

        // Close the mat file
        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }
        cout<<"--------finish writing--------"<<endl;
    }


    /*
     * create a 3d matrix
     */
    float*** create3dMatrix_float(int h, int w, int l){
        float*** data;
        data = new float** [h];
        for (int i = 0; i < h; ++i) {
            data[i] = new float* [w];
        }
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                data[i][j] = new float [l];
            }
        }
//        for (int i = 0; i < h; ++i) {
//            for (int j = 0; j < w; ++j) {
//                for (int k = 0; k < l; ++k) {
//                    data[i][j][k] = 0;
//                }
//            }
//        }
        return data;
    }// create3dMatrix_float

    int*** create3dMatrix_int(int h, int w, int l){
        int*** data;
        data = new int** [h];
        for (int i = 0; i < h; ++i) {
            data[i] = new int* [w];
        }
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                data[i][j] = new int [l];
            }
        }
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                for (int k = 0; k < l; ++k) {
                    data[i][j][k] = 0;
                }
            }
        }
        return data;
    }// create3dMatrix_int

    int**** create4dMatrix_int(int h, int w, int l, int t){
        int**** data;
        data = new int*** [h];
        for (int i = 0; i < h; ++i) {
            data[i] = new int** [w];
        }
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                data[i][j] = new int* [l];
            }
        }
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                for (int k = 0; k < l; ++k) {
                    data[i][j][k] = new int [t];
                }
            }
        }
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                for (int k = 0; k < l; ++k) {
                    for (int tt = 0; tt < t; ++tt) {
                        data[i][j][k][tt] = 0;
                    }
                }
            }
        }
        return data;
    }// create4dMatrix()


    bool*** createEvtSpatialMask(int H, int W, int L){
        bool*** evtSpatialMask;
        evtSpatialMask = new bool** [H];
        //#pragma omp parallel for
        for (int i = 0; i < H; ++i) {
            evtSpatialMask[i] = new bool* [W];
        }
        //#pragma omp parallel for collapse(2)
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                evtSpatialMask[i][j] = new bool [L];
            }
        }
        //#pragma omp parallel for collapse(3)
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    evtSpatialMask[i][j][k] = true;
                }
            }
        }
        return evtSpatialMask;
    }


    void release3dMatrix(float***& data, int h, int w){
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                delete[] data[i][j];
            }
            delete[] data[i];
        }
        delete[] data;
        data = nullptr;
    } // release3dMatrix


    void release3dMatrix_bool(bool***& data, int h, int w){
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                delete[] data[i][j];
            }
            delete[] data[i];
        }
        delete[] data;
        data = nullptr;
    } // release3dMatrix_bool

    
    void release3dMatrix_int(int***& data, int h, int w){
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w; ++j) {
                delete[] data[i][j];
            }
            delete[] data[i];
        }
        delete[] data;
        data = nullptr;
    } // release3dMatrix_int

    int sub2ind(int i, int j, int k, int h, int w){
        return i + j*h + k*h*w;
    }

    int sub2ind(int i, int j, int k,int t, int h, int w, int l){
        return i + j*h + k*h*w + t*h*w*l;
    }

    Point_struct ind2sub(int ind, int h, int w){ // column first
        Point_struct ans;
        ans.k = ind / (h*w);
        ind -= ans.k * h * w;
        ans.j = ind / h;
        ans.i = ind % h;
        return ans;
    }

    Point_struct ind2sub(int ind, int h, int w, int l){ // column first
        Point_struct ans;
        ans.t = ind / (h*w*l);
        ind -= ans.t * h * w * l;
        ans.k = ind / (h*w);
        ind -= ans.k*h*w;
        ans.j = ind / h;
        ans.i = ind % h;
        return ans;
    }


    void save_vector(const vector<double>& vec, const string& file_path) {
        ofstream output_file(file_path, ios::binary);
        if (!output_file.is_open()) {
            throw runtime_error("Unable to open file for writing: " + file_path);
        }
        
        size_t size = vec.size();
        output_file.write(reinterpret_cast<const char*>(&size), sizeof(size));
        
        output_file.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(double));
        output_file.close();
        cout<<"vector saved"<< endl;
    }


    
    vector<double> load_vector(const string& file_path) {
        ifstream input_file(file_path, ios::binary);
        if (!input_file.is_open()) {
            throw runtime_error("Unable to open file for reading: " + file_path);
        }
        
        size_t size;
        input_file.read(reinterpret_cast<char*>(&size), sizeof(size));
        
        vector<double> vec(size);
        input_file.read(reinterpret_cast<char*>(vec.data()), size * sizeof(double));
        input_file.close();
        cout<<"vector loaded"<< endl;
        return vec;
    }
    




}// namespace