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
            for (int k = 0; k < L; ++k) {
                frame[t].emplace_back(H- 2*bdCrop,W- 2*bdCrop,CV_32F);
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
        std::cout<<"after cropping: "<< std::endl;
        std::cout<<"height of image:"<< H << std::endl;
        std::cout<<"width of image:"<< W << std::endl;
        std::cout<<"length of image:"<< L << std::endl;
        std::cout<<"time frames of image:"<< T << std::endl;

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
        std::cout<<"--------data loaded--------"<<std::endl;

        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        frame[t][k].at<float>(i,j) = (frame[t][k].at<float>(i,j) - static_cast<float>(AQuA::opts.minValueDat1)) / normalizedParameter;
//                    std::cout<< frame[t][k].at<float>(i,j)<< "  ";// display pixel value
                    }//for(i)
                }//for(j)
            }//for(k)
        }//for(t)
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
//
        return frame;
    }//loadData()


    /*
     * create a 3d matrix
     */
    float*** create3dMatrix(int h, int w, int l){
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
//                for (int k = 0; k < L; ++k) {
//                    data[i][j][k] = 0;
//                }
//            }
//        }
        return data;
    }// create3dMatrix()


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
        opts.BitDepth;
        opts.regMaskGap = 5;
        opts.singleChannel = true;
        opts.registrateCorrect = 0;
        opts.bleachCorrect = 0;
        opts.medSmo = 1;
        opts.cut = 200;

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
//    /*
//     * load image data to a std::vector<std::vector<cv::Mat>>structure
//     * access pixel value by frame[t][k].at<float>(i,j)
//     */
//    std::vector<std::vector<cv::Mat>> loadData() {
//        T = 1;
//        std::vector<std::vector<cv::Mat>> frame(T);
//        std::string pre_name = "C:/Users/Kevin Qiao/Desktop/3D_data/3D_dataFrame ";
//        std::string name_ext = ".tif";
//        double min,max=0;
//        double mmin=255, mmax=0;
//        int bdCrop = opts.regMaskGap;
//        int BitDepth = -1;
//        float normalizedParameter;
//
//        std::cout<<"--------reading data--------"<<std::endl;
//        for (int t = 0; t < T; ++t) {
//            cv::imreadmulti(pre_name + std::to_string(t+1) + name_ext, frame[t],  cv::IMREAD_GRAYSCALE); // uchar 8-bit unsigned int
//            if (t == 0) {
//                H = frame[0][0].rows;
//                W = frame[0][0].cols;
//                L = static_cast<int>(frame[0].size());
//            }//(t=0) output size of each dimension
//            for (int k = 0; k < L; ++k) {
//                crop(frame[t][k],bdCrop);
//                if (t==0 & k==0) {
//                    H = frame[0][0].rows;
//                    W = frame[0][0].cols;
//                    std::cout<<"after cropping:"<< std::endl;
//                    std::cout<<"height of image:"<< H << std::endl;
//                    std::cout<<"width of image:"<< W << std::endl;
//                    std::cout<<"length of image:"<< L << std::endl;
//                }//(t=0) output size of each dimension
//                cv::minMaxLoc(frame[t][k], &min, &max);
//                if (min < mmin){
//                    mmin = min;
//                }
//                if(max > mmax){
//                    mmax = max;
//                }
////            std::cout<< min<< " "<<max<< std::endl;
////            cv::imwritemulti("C:/Users/Kevin Qiao/Desktop/1.tif",frame[t]); //save image
////            std::cout << "slice " << k + 1 << " of " << L << " frame " << t + 1 << " of " << T << std::endl;
////            cv::imshow("length", frame[t][k]);
////            cv::waitKey();
//            }//for(k)
//        }//for(t)
//        std::cout<<"--------data loaded--------"<<std::endl;
//        AQuA::opts.maxValueDat1 = mmax;
//        AQuA::opts.minValueDat1 = mmin;
//        normalizedParameter = static_cast<float>(mmax -mmin);
//        for (int t = 0; t < T; ++t) {
//            for (int k = 0; k < L; ++k) {
//                frame[t][k].convertTo(frame[t][k], CV_32F);
//                for (int i = 0; i < H; ++i) {
//                    for (int j = 0; j < W; ++j) {
//                        frame[t][k].at<float>(i,j) -= static_cast<float>(AQuA::opts.minValueDat1);
////                    if (frame[t][k].at<float>(i,j)<0) {
////                        frame[t][k].at<float>(i, j) = 0;
////                    }
//                        frame[t][k].at<float>(i,j) /= normalizedParameter;
////                    std::cout<< frame[t][k].at<float>(i,j)<< "  ";// display pixel value
//                    }//for(i)
//                }//for(j)
////            std::cout << "slice " << k + 1 << " of " << L << " frame " << t + 1 << " of " << T << std::endl;
////            cv::imshow("length", frame[t][k]);
////            cv::waitKey();
////            std::cout<<std::endl;
//            }//for(k)
////        std::cout<<std::endl;
//        }//for(t)
//        AQuA::opts.sz[0] = H;
//        AQuA::opts.sz[1] = W;
//        AQuA::opts.sz[2] = L;
//        AQuA::opts.sz[3] = T;
//        AQuA::opts.BitDepth = BitDepth;
//        return frame;
//    }//loadData()

//    /*
//     * create a 4d matrix, which size is T,L,H,W
//     * initialize the values with 0
//     */
//    DATA_TYPE**** create4dMatrix(){
//
//        DATA_TYPE**** data;
//        data = new DATA_TYPE*** [H];
//        for (int i = 0; i < H; ++i) {
//            data[i] = new DATA_TYPE** [W];
//
//        }
//        for (int i = 0; i < H; ++i) {
//            for (int j = 0; j < W; ++j) {
//                data[i][j] = new DATA_TYPE* [L];
//            }
//        }
//        for (int i = 0; i < H; ++i) {
//            for (int j = 0; j < W; ++j) {
//                for (int k = 0; k < L; ++k) {
//                    data[i][j][k] = new DATA_TYPE [T];
//                }
//            }
//        }
//
////        for (int i = 0; i < H; ++i) {
////            for (int j = 0; j < W; ++j) {
////                for (int k = 0; k < L; ++k) {
////                    for (int t = 0; t < T; ++t) {
////                        data[i][j][k][t] = 0;
////                    }
////                }
////            }
////        }
//
//
//        return data;
//    }// getData()
//
//

//
//
//    /*
//     * create a 3d matrix which size is H_ext,W_ext,L_ext
//     * initialize the values with 0
//     */
//    double *** create3dMatrix_ext_double(){
//        double*** data;
//        data = new double ** [H_ext];
//        for (int i = 0; i < H_ext; ++i) {
//            data[i] = new double * [W_ext];
//        }
//        for (int i = 0; i < H_ext; ++i) {
//            for (int j = 0; j < W_ext; ++j) {
//                data[i][j] = new double [L_ext];
//            }
//        }
////        for (int i = 0; i < H_ext; ++i) {
////            for (int j = 0; j < W_ext; ++j) {
////                for (int k = 0; k < L_ext; ++k) {
////                    data[i][j][k] = 0;
////                }
////            }
////        }
//        return data;
//    }// create3dMatrix_ext()
//
//
//    float *** create3dMatrix_ext_float(){
//        float *** data;
//        data = new float ** [H_ext];
//        for (int i = 0; i < H_ext; ++i) {
//            data[i] = new float * [W_ext];
//        }
//        for (int i = 0; i < H_ext; ++i) {
//            for (int j = 0; j < W_ext; ++j) {
//                data[i][j] = new float [L_ext];
//            }
//        }
//        for (int i = 0; i < H_ext; ++i) {
//            for (int j = 0; j < W_ext; ++j) {
//                for (int k = 0; k < L_ext; ++k) {
//                    data[i][j][k] = 0;
//                }
//            }
//        }
//        return data;
//    }// create3dMatrix_ext()_float()
//
//
//
//    void releaseData(int* data){
//        delete[] data;
//        data = NULL;
//    }// releaseData
//
//    void releaseData(DATA_TYPE* data){
//        delete[] data;
//        data = NULL;
//    }// releaseData
//
//    void releaseData(DATA_TYPE**** data, int I, int J, int K){
//        for (int i = 0; i < I; ++i) {
//            for (int j = 0; j < J; ++j) {
//                for (int k = 0; k < K; ++k) {
//                    delete[] data[i][j][k];
//                }
//                delete[] data[i][j];
//            }
//            delete[] data[i];
//        }
//        delete[] data;
//        data = NULL;
//    }// releaseData


//    void releaseData(double*** data, int I, int J){
//        for (int i = 0; i < I; ++i) {
//            for (int j = 0; j < J; ++j) {
//                delete[] data[i][j];
//            }
//            delete[] data[i];
//        }
//        delete[] data;
//        data = NULL;
//    }// releaseData






}// namespace