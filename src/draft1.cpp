//#include <opencv2/opencv.hpp>
//#include "draft.h"
#include "data/data.h"
#include <chrono>

namespace AQuA{

    float medianFunc(float* array, int size){
        float median;
        std::sort(array, array+size);
        if(size % 2 == 0){
            median = ( array[size/2] + array[size/2 - 1] ) / 2;
        }
        else{
            median = array[size/ 2];
        }
        std::cout<<"median value of reference frame: "<< median<<std::endl;
        return median;
    }//medianFunc()


    float*** flip3dMatrix(float*** ref, float*** b_flip){
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    b_flip[H - 1 - i][W - 1 -j ][L - 1 - k] = ref[i][j][k];
                }// for(k)
            }// for(j)
        }// for(i)
        return b_flip;
    }// flip_3d()
    
    
    float*** calCC(float*** a, float*** b, float*** a_add, float*** b_add, float*** b_flip){
        b_flip = flip3dMatrix(b,b_flip); // flip(flip(flip(b,1),2),3);
        for(int i=0;i<H;++i){
            for(int j=0;j<W;++j){
                for (int k = 0; k < L; ++k) {
                    a_add[i][j][k] = a[i][j][k];
                    b_add[i][j][k] = b_flip[i][j][k];
                } //for(k)
            } //for(j)
        }// for(i)
        return dft(a_add,b_add);
    }// calCC


    void regCrossCorrelation(std::vector<std::vector<cv::Mat>>& data1){
        float mean_sum, median;
        int refer_start=0, refer_end =9;
        float*** ref = create3dMatrix(H,W,L); //remember to release with release3dMatrix(), the following matrix as well
        float* array_1d = new float [H*W*L]; //remember to release with delete[]
        float*** moving = create3dMatrix(H,W,L);
        float*** a_add = create3dMatrix(H_ext,W_ext,L_ext);
        float*** b_add = create3dMatrix(H_ext,W_ext,L_ext);
        float*** b_flip = create3dMatrix(H,W,L);
        float*** matrix = create3dMatrix(H,W,L);
        float* matrix_col = new float [H_ext*W_ext*L_ext];

        /*
         * calculate mean value in each time frame and store in a 3d matrix ref[][][]
         */
        for (int i=0, m=0; i<H; ++i) {
            for (int j=0; j<W; ++j) {
                for(int k=0; k<L;++k){
                    mean_sum = 0;
                    for (int t=refer_start;t<=refer_end;++t) {
                        mean_sum += data1[t][k].at<float>(i,j);
                    }// for(t)
                    ref[i][j][k] = mean_sum / static_cast<float>(refer_end - refer_start + 1);
                    array_1d[m++] = ref[i][j][k]; //convert ref[][][] to one dimension, 'm' is just index for the array
                }// for(k)
            }// for(j)
        }// for(i)

        median = medianFunc(array_1d, H*W*L);

        /*
         * align bright part. Remove median is like remove background
         */
        for(int i=0;i<H;++i){
            for(int j=0;j<W;++j){
                for(int k=0;k<L;++k){
                    ref[i][j][k] -= median;  // ref = ref - median(ref(:));
                }// for(k)
            }// for(j)
        }// for(i)

        /*
         * concise cross correlation
         */
        for(int t=0; t<T; ++t) {
            for (int i = 0, m = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    for (int k = 0; k < L; ++k) {
                        moving[i][j][k] = data1[t][k].at<float>(i, j); //initialize moving
                        array_1d[m++] = moving[i][j][k]; //convert moving[][][] to one dimension
                    }// for(k)
                } // for(j)
            }// for(i)
            median = medianFunc(array_1d, H * W * L);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    for (int k = 0; k < L; ++k) {
                        moving[i][j][k] -= median;  // moving = moving - median(moving(:));
                        matrix = calCC(moving, ref, a_add, b_add, b_flip); // matrix = calCC(moving,ref);
                    }// for(k)
                } //for(j)
            }// for(i)


            release3dMatrix(ref, H, W);
            delete[] array_1d;
        }//for(t)
        
    }//regCrossCorrelation()



}//namespace


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



