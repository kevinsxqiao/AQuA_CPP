//
// Created by Kevin Qiao on 3/28/2023.
//

#include "data.h"


namespace AQuA{

    /*
     * define data structure
     */
    struct rawDataSize_struct rawDataSize;
    struct preSetting_struct preSetting;
    struct opts_struct opts;

    std::vector<std::vector<cv::Mat>> loadData() {
        MATFile *pmatFile;
        mxArray *pMxArray;
        double *pdata;
        int bdCrop = AQuA::opts.regMaskGap;
        double min,max=0;
        double mmin=255, mmax=0;
        int BitDepth = -1;
        float normalizedParameter;

        std::cout<< "--------loading data--------"<<std::endl;
//        std::cout<<'\r'<< 0 << "% "<<std::flush;
//        const char *filename = "C:/Users/Kevin Qiao/Desktop/AQuA_data/Test_global_local_3D.mat";
        pmatFile = matOpen(opts.fileName1, "r");
        if (pmatFile == nullptr) {
            std::cout<< "--------error opening file--------"<<std::endl;
            std::exit(-1);
        }

        pMxArray = matGetVariable(pmatFile, "synthetic");
        if (pMxArray == nullptr) {
            std::cout<< "--------error reading variable from file--------"<<std::endl;
            std::exit(-1);
        }

        pdata = mxGetPr(pMxArray);
        if (pdata == nullptr) {
            std::cout<< "--------error reading data from variable-------"<<std::endl;
            std::exit(-1);
        }

        const mwSize *dims = mxGetDimensions(pMxArray);
        H = dims[0];
        W = dims[1];
        L = dims[2];
        T = dims[3];
//        T = 10;
//        std::cout<<"original size: "<< std::endl;
//        std::cout<<"height of image:"<< H << std::endl;
//        std::cout<<"width of image:"<< W << std::endl;
//        std::cout<<"length of image:"<< L << std::endl;
//        std::cout<<"time frames of image:"<< T << std::endl;

        std::vector<std::vector<cv::Mat>> frame(T);
        for (int t = 0; t < T; ++t) {
//            std::cout<<'\r'<< 100*t/T/2 <<"% "<<std::flush;
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
//        std::cout<<"minValue: "<< opts.minValueDat1 << "  maxValue: "<<opts.maxValueDat1<<std::endl;
        normalizedParameter = static_cast<float>(mmax -mmin);

//#pragma omp parallel for collapse(4)
        for (int t = 0; t < T; ++t) {
//            std::cout<<'\r'<< 50+ (100*t/T/2) <<"%"<<std::flush;
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        frame[t][k].at<float>(i,j) = (frame[t][k].at<float>(i,j) - static_cast<float>(AQuA::opts.minValueDat1)) / normalizedParameter;
//                    std::cout<< frame[t][k].at<float>(i,j)<< "  ";// display pixel value
                    }//for(i)
                }//for(j)
            }//for(k)
        }//for(t)
//        std::cout << "\r" << 100 << "%" << std::flush<< std::endl;
        AQuA::opts.sz[0] = H;
        AQuA::opts.sz[1] = W;
        AQuA::opts.sz[2] = L;
        AQuA::opts.sz[3] = T;
        AQuA::opts.BitDepth = BitDepth;

//        std::cout<<"data: "<<std::endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                std::cout<< frame[0][0].at<float>(i,j)<< "  ";
//            }
//            std::cout<<std::endl;
//        }
    std::cout<<"after cropping: "<< std::endl;
    std::cout<<"height of image:"<< H << std::endl;
    std::cout<<"width of image:"<< W << std::endl;
    std::cout<<"length of image:"<< L << std::endl;
    std::cout<<"time frames of image:"<< T << std::endl;
    std::cout<<"--------data loaded--------"<<std::endl;
    opts.data1_org = frame;

//    for (int t = 0; t < T; ++t) {
//        std::vector<cv::Mat> dat(L);
//        for (int k = 0; k < L; ++k) {
//            dat[t].push_back(frame[t][k].clone());
//        }
//        opts.data1_org.push_back(dat);
//    }

    return frame;
    }//loadData()


    std::vector<std::vector<cv::Mat>> load4DData_clean(const char* fileName, const char* varName) {
        MATFile *pmatFile;
        mxArray *pMxArray;
        double *pdata;

        std::cout<< "--------loading data--------"<<std::endl;
        pmatFile = matOpen(fileName, "r");
        if (pmatFile == nullptr) {
            std::cout<< "--------error opening file--------"<<std::endl;
            std::exit(-1);
        }

        pMxArray = matGetVariable(pmatFile, varName);
//        pMxArray = matGetVariable(pmatFile, "dF");
        if (pMxArray == nullptr) {
            std::cout<< "--------error reading variable from file--------"<<std::endl;
            std::exit(-1);
        }

        pdata = mxGetPr(pMxArray);
        if (pdata == nullptr) {
            std::cout<< "--------error reading data from variable-------"<<std::endl;
            std::exit(-1);
        }

        const mwSize *dims = mxGetDimensions(pMxArray);
        H = dims[0];
        W = dims[1];
        L = dims[2];
        T = dims[3];
        std::cout<<"original size: "<< std::endl;
        std::cout<<"height of image:"<< H << std::endl;
        std::cout<<"width of image:"<< W << std::endl;
        std::cout<<"length of image:"<< L << std::endl;
        std::cout<<"time frames of image:"<< T << std::endl;

        std::vector<std::vector<cv::Mat>> frame(T,std::vector<cv::Mat>(L));
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                frame[t][k] = cv::Mat(H,W,CV_32F);
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j){
                        frame[t][k].at<float>(i,j) = static_cast<float>(pdata[sub2ind(i,j,k,t,H,W,L)]);
//                        std::cout<<t<<" "<<k<<" "<<i<<" "<<j<<" "<<std::endl;
                    }//for(j)
                }//for(i)
            }//for(k)
        }//for(t)


        //release MAT pointer
        if (pMxArray != nullptr) {
            mxDestroyArray(pMxArray);
        }

        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }

        std::cout<<"height of image:"<< H << std::endl;
        std::cout<<"width of image:"<< W << std::endl;
        std::cout<<"length of image:"<< L << std::endl;
        std::cout<<"time frames of image:"<< T << std::endl;
        std::cout<<"--------data loaded--------"<<std::endl;
        return frame;
    }//loadData()


    mxArray* cvDataToMxArray(const std::vector<std::vector<cv::Mat>>& data) {
        // Calculate the size of the 4D matrix
        mwSize dims[4] = {static_cast<mwSize>(data[0][0].rows), static_cast<mwSize>(data[0][0].cols), static_cast<mwSize>(data[0].size()), static_cast<mwSize>(data.size())};

        // Create a 4D mxArray
        mxArray* pMxArray = mxCreateNumericArray(4, dims, mxSINGLE_CLASS, mxREAL);

        // Copy data from your vector to the mxArray
        float* ptr = reinterpret_cast<float*>(mxGetData(pMxArray));
        for (int t = 0; t < data.size(); ++t) {
            for (int k = 0; k < data[t].size(); ++k) {
                const cv::Mat& mat = data[t][k];
                std::memcpy(ptr, mat.data, mat.rows * mat.cols * sizeof(float));
                ptr += mat.rows * mat.cols;
            }
        }

        return pMxArray;
    }


    mxArray* cvDataToMxArray(const std::vector<cv::Mat>& data) {
        // Calculate the size of the 3D matrix
        mwSize dims[3] = {static_cast<mwSize>(data[0].rows), static_cast<mwSize>(data[0].cols), static_cast<mwSize>(data.size())};

        // Create a 3D mxArray
        mxArray* pMxArray = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);

        // Copy data from your vector to the mxArray
        float* ptr = reinterpret_cast<float*>(mxGetData(pMxArray));
            for (int k = 0; k < data.size(); ++k) {
                const cv::Mat& mat = data[k];
                std::memcpy(ptr, mat.data, mat.rows * mat.cols * sizeof(float));
                ptr += mat.rows * mat.cols;
            }
        return pMxArray;
    }


    void writeDataToMatFile(std::vector<std::vector<cv::Mat>>& data, const std::string& filename) {
        std::cout<<"--------start writing--------"<<std::endl;
        MATFile *pmatFile;

        // Open the mat file
        pmatFile = matOpen(filename.c_str(), "w");
        if (pmatFile == nullptr) {
            std::cout << "--------error opening file--------" << std::endl;
            std::exit(-1);
        }

        // Convert your data to a mxArray
        mxArray* pMxArray = cvDataToMxArray(data);

        // Write the variable to the mat file
        if (matPutVariable(pmatFile, "myVar", pMxArray) != 0) {
            std::cout << "--------error writing variable to file--------" << std::endl;
            std::exit(-1);
        }

        // Free the mxArray
        mxDestroyArray(pMxArray);

        // Close the mat file
        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }
        std::cout<<"--------finish writing--------"<<std::endl;
    }


    void writeDataToMatFile(std::vector<cv::Mat>& data, const std::string& filename) {
        std::cout<<"--------start writing--------"<<std::endl;
        MATFile *pmatFile;

        // Open the mat file
        pmatFile = matOpen(filename.c_str(), "w");
        if (pmatFile == nullptr) {
            std::cout << "--------error opening file--------" << std::endl;
            std::exit(-1);
        }

        // Convert your data to a mxArray
        mxArray* pMxArray = cvDataToMxArray(data);

        // Write the variable to the mat file
        if (matPutVariable(pmatFile, "myVar", pMxArray) != 0) {
            std::cout << "--------error writing variable to file--------" << std::endl;
            std::exit(-1);
        }

        // Free the mxArray
        mxDestroyArray(pMxArray);

        // Close the mat file
        if (pmatFile != nullptr) {
            matClose(pmatFile);
        }
        std::cout<<"--------finish writing--------"<<std::endl;
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


    bool*** createEvtSpatialMask(){
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


//    judge if registration and bleach have been executed; ---- true = no ; false = both executed
    bool isDefault() {
        return (preSetting.registrateCorrect == preSetting.registrateCorrect_default
                && preSetting.bleachCorrect == preSetting.bleachCorrect_default);
    }// isDefault()


    void preSettingInit(){
        preSetting.registrateCorrect = preSetting.registrateCorrect_default;
        preSetting.bleachCorrect = preSetting.bleachCorrect_default;
        std::cout<< "--------preSetting initialized--------"<<std::endl;
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

        std::cout<< "--------opts initialized--------"<<std::endl;
    }// optsInit()


    void rawDataSizeInit(){
        rawDataSize={0};
        std::cout<< "--------rawDataSize initialized--------"<<std::endl;
    }


    void Init(){
        preSettingInit();
        optsInit();
        rawDataSizeInit();
        std::cout<<std::endl;
    }

//    /*
//     * crop the 2d image from each direction by 'bdCrop'
//     */
//    void crop(cv::Mat& image, int bdCrop){
////    cv::Rect roi = cv::Rect(bdCrop,bdCrop,W-2*bdCrop,H-2*bdCrop);
////    image.adjustROI(-roi.y, image.rows - (roi.y+roi.height), -roi.x, image.cols- (roi.x+roi.width));
//        image.adjustROI(-bdCrop, -bdCrop, -bdCrop, -bdCrop);
//    }
//
//




}// namespace