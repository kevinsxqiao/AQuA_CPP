//
// Created by Kevin Qiao on 2/17/2023.
//

#include "regCrossCorrelation.h"

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
//        std::cout<<"median value of reference frame: "<< median<<std::endl;
        return median;
    }//medianFunc()

    /*
     *   allocate memory  --use recommend fftwf_malloc which allocating memory on the heap
     *   fixed-size array is allocated on the stack, using large fixed-size arrays may cause stack overflow issues or exceed the available stack size.
     */
    float*** dft(float*** a_add, float*** b_add){
        float* input_a = static_cast<float *>(fftwf_malloc(H_ext * W_ext * L_ext * sizeof(float)));
        float* input_b = static_cast<float *>(fftwf_malloc(H_ext * W_ext * L_ext * sizeof(float)));
        fftwf_complex* output_a = static_cast<fftwf_complex *>(fftwf_malloc(H_ext * W_ext * (L_ext/2+1) * sizeof(fftwf_complex)));
        fftwf_complex* output_b = static_cast<fftwf_complex *>(fftwf_malloc(H_ext * W_ext * (L_ext/2+1) * sizeof(fftwf_complex)));
        fftwf_complex* output_dotProduct = static_cast<fftwf_complex *>(fftwf_malloc(H_ext * W_ext * (L_ext/2+1) * sizeof(fftwf_complex)));
        float* result = static_cast<float *>(fftwf_malloc(H_ext * W_ext * L_ext * sizeof(float)));
        float*** c = create3dMatrix(H_ext,W_ext,L_ext);

        // flatten the 3d image data and store in input_a[]
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                for (int k = 0; k < L_ext; ++k) {
                    int index = i * W_ext * L_ext + j * L_ext + k;
                    input_a[index] = a_add[i][j][k];
                    input_b[index] = b_add[i][j][k];
//                    std::cout<<input_a[i * W_ext * L_ext + j * L_ext + k]<<" ";
                }
            }
        }

//        std::cout<<"input_a"<<std::endl;
//            for (int i = 0; i < 10; ++i) {
//                for (int j = 0; j < 10; ++j) {
//                    std::cout<<input_a[i * W_ext * L_ext + j * L_ext ]<<" ";
//                }
//                std::cout<<std::endl;
//            }
//            std::cout<<std::endl;

        // create fftw plan
        fftwf_plan fft_plan_a = fftwf_plan_dft_r2c_3d(H_ext, W_ext, L_ext, input_a, output_a, FFTW_ESTIMATE);
        fftwf_plan fft_plan_b = fftwf_plan_dft_r2c_3d(H_ext, W_ext, L_ext, input_b, output_b, FFTW_ESTIMATE);
        fftwf_plan ifft_plan = fftwf_plan_dft_c2r_3d(H_ext, W_ext, L_ext, output_dotProduct, result,  FFTW_ESTIMATE);

        fftwf_execute(fft_plan_a);
        fftwf_execute(fft_plan_b);

////        display the matrix after fft_plan_a
//        std::cout<<"output_a"<<std::endl;
//        for (int k = 0; k < 1; ++k) {
//            for (int i = 0; i < 10; ++i) {
//                for (int j = 0; j < 10; ++j) {
//                    int index = i * W_ext * (L_ext/2 + 1) + j * (L_ext/2 + 1) + k;
//                    std::cout<< output_a[index][0]<< "+" << output_a[index][1]<< "i"<< "      ";
//                }
//                std::cout<< std::endl;
//            }
//            std::cout<< std::endl;
//        }

        // dot product
        for (int index = 0; index<H_ext * W_ext * (L_ext/2+1); ++index) {
            output_dotProduct[index][0] = output_a[index][0] * output_b[index][0] - output_a[index][1] * output_b[index][1];
            output_dotProduct[index][1] = output_a[index][0] * output_b[index][1] + output_a[index][1] * output_b[index][0];
//            std::cout<< output_dotProduct[index][0]<< "+" << output_dotProduct[index][1]<< "i"<< "      ";
        }

//        std::cout<<"output_dotProduct"<<std::endl;
//        for (int k = 0; k < 1; ++k) {
//            for (int i = 0; i < 10; ++i) {
//                for (int j = 0; j < 10; ++j) {
//                    int index = i * W_ext * (L_ext/2 + 1) + j * (L_ext/2 + 1) + k;
//                    std::cout<< output_dotProduct[index][0]<< "+" << output_dotProduct[index][1]<< "i"<< "      ";
//                }
//                std::cout<< std::endl;
//            }
//            std::cout<< std::endl;
//        }
//
//        fftwf_execute(ifft_plan);



        // normalize the IFFT output, since ifftn() will automatically normalize in matlab
        int total = H_ext*W_ext*L_ext;
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                for (int k = 0; k < L_ext; ++k) {
                    int index = i * W_ext * L_ext + j * L_ext + k;
                    c[i][j][k] = result[index]/static_cast<float>(total) - 9.45;//difference between library
//                    std::cout<< c[i][j][k]<< " ";
                }
            }
        }
//        std::cout<<" c=ifftn(fftn(a_add).*fftn(b_add));"<<std::endl;
//        for (int i = 0; i < 10; ++i) {
//            for (int j = 0; j < 10; ++j) {
//                std::cout<<c[i][j][0]<< " ";
//            }
//            std::cout<<std::endl;
//        }

        // release memory
        fftwf_destroy_plan(fft_plan_a);
        fftwf_destroy_plan(fft_plan_b);
        fftwf_destroy_plan(ifft_plan);
        fftwf_free(input_a);
        fftwf_free(input_b);
        fftwf_free(output_a);
        fftwf_free(output_b);
        fftwf_free(result);
        fftwf_free(output_dotProduct);

        return c;
    }//dft()


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

//        for (int i = 0; i < 10; ++i) {
//            for (int j = 0; j < 10; ++j) {
//                std::cout<<a_add[i][j][0]<< " ";
//            }
//            std::cout<<std::endl;
//        }
//        std::cout<<std::endl;
//        for (int i = 0; i < 10; ++i) {
//            for (int j = 0; j < 10; ++j) {
//                std::cout<<b_add[i][j][0]<< " ";
//            }
//            std::cout<<std::endl;
//        }
//        std::cout<<std::endl;
        return dft(a_add,b_add);
    }// calCC


    std::vector<std::vector<cv::Mat>> regCrossCorrelation(std::vector<std::vector<cv::Mat>>& data1){
        float mean_sum, median;
        int refer_start=0, refer_end =9;
        float*** ref = create3dMatrix(H,W,L); //remember to release with release3dMatrix(), the following matrix as well
        float* array_1d = new float [H*W*L]; //remember to release with delete[]
        float*** moving = create3dMatrix(H,W,L);
        float*** a_add = create3dMatrix(H_ext,W_ext,L_ext);
        float*** b_add = create3dMatrix(H_ext,W_ext,L_ext);
        float*** b_flip = create3dMatrix(H,W,L);
        float*** matrix = create3dMatrix(H_ext,W_ext,L_ext);
//        float* matrix_1d = new float [H_ext*W_ext*L_ext];
//        int id,r1;
        float maxElement = -100000;
        int wShift, hShift, lShift;
        int xs0, xe0, xs1, xe1, ys0, ye0, ys1, ye1, zs0, ze0, zs1, ze1;
        int* x_translation = new int [T];
        int* y_translation = new int [T];
        int* z_translation = new int [T];

        std::cout<<"--------start regCrossCorrelation--------"<<std::endl;
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
//                    std::cout<<ref[i][j][k]<< " ";
                    array_1d[m++] = ref[i][j][k]; //flatten ref[], 'm' is just index for the array
                }// for(k)
            }// for(j)
        }// for(i)

        median = medianFunc(array_1d, H*W*L);

        std::cout<<"--------start align bright part--------"<<std::endl;
        /*
         * align bright part. Remove median is like remove background
         */
        for(int i=0, m=0;i<H;++i){
            for(int j=0;j<W;++j){
                for(int k=0;k<L;++k){
                    ref[i][j][k] -= median;  // ref = ref - median(ref(:));
                }// for(k)
            }// for(j)
        }// for(i)

//        for (int i = 0; i < 10; ++i) {
//            for (int j = 0; j < 10; ++j) {
//                std::cout<<ref[i][j][0]<< " ";
//            }
//            std::cout<<std::endl;
//        }

        std::cout<<"--------start cross correlation--------"<<std::endl;
        /*
         * concise cross correlation
         */
        for(int t=0; t<T; ++t) {
            for (int i = 0, m = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    for (int k = 0; k < L; ++k) {
                        moving[i][j][k] = data1[t][k].at<float>(i, j); //initialize moving
                        array_1d[m++] = moving[i][j][k]; //flatten moving[]
                    }// for(k)
                } // for(j)
            }// for(i)
            median = medianFunc(array_1d, H * W * L);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    for (int k = 0; k < L; ++k) {
                        moving[i][j][k] -= median;  // moving = moving - median(moving(:));
                    }// for(k)
                } //for(j)
            }// for(i)
//            for (int i = 0; i < 10; ++i) {
//                for (int j = 0; j < 10; ++j) {
//                    std::cout<<moving[i][j][0]<< " ";
//                }
//                std::cout<<std::endl;
//            }

//            std::cout<<"--------start calCC--------"<<std::endl;
            matrix = calCC(moving, ref, a_add, b_add, b_flip); // matrix = calCC(moving,ref); run the first time is time-consuming, try 'wisdom' in fftw
//            std::cout<<"--------end calCC--------"<<std::endl;

//            /*
//             * flatten matrix[] to find the position of the maximum element
//             */
//            for(int i=0, m=0;i<H_ext;++i){
//                for(int j=0;j<W_ext;++j){
//                    for(int k=0;k<L_ext;++k){
//                        matrix_1d[m++] = matrix[i][j][k];
//                    }// for(k)
//                } //for(j)
//            }// for(i)
//            id = std::max_element(matrix_1d, matrix_1d+(H_ext*W_ext*L_ext)) - matrix_1d;  //[~,id] = max(matrix(:));
            for(int i=0;i<H_ext;++i){
                for(int j=0;j<W_ext;++j){
                    for(int k=0;k<L_ext;++k){
//                        std::cout<<matrix[i][j][k]<<" ";
                        if (matrix[i][j][k] > maxElement){
                            maxElement = matrix[i][j][k];
                            hShift = i+1;
                            wShift = j+1;
                            lShift = k+1;
//                            std::cout<<"hShift: "<<hShift<<"wShift: "<<wShift<<"lShift: "<<lShift<<std::endl;
                        }
                    }// for(k)
                } //for(j)
            }// for(i)
            // [hShift,wShift,lShift] = ind2sub(size(matrix),id);
//            hShift = id/(W_ext*L_ext);
//            r1 = id%(W_ext*L_ext);
//            wShift = r1 / L_ext;
//            lShift = r1 % L_ext;
            x_translation[t] = H - hShift;
            y_translation[t] = W - wShift;
            z_translation[t] = L - lShift;
        }//for(t)

        std::cout<<"--------start translation--------"<<std::endl;
        for(int t=0;t<T;++t) {
            if (x_translation[t] >= 0) {//towards right
                xs0 = 0;
                xe0 = H - x_translation[t] - 1;
                xs1 = x_translation[t];
                xe1 = H - 1;
            } else {//towards left
                xs0 = -x_translation[t];
                xe0 = H - 1;
                xs1 = 0;
                xe1 = H + x_translation[t] - 1;
            }

            if (y_translation[t] >= 0) {
                ys0 = 0;
                ye0 = W - y_translation[t] - 1;
                ys1 = y_translation[t];
                ye1 = W - 1;
            } else {
                ys0 = -y_translation[t];
                ye0 = W - 1;
                ys1 = 0;
                ye1 = W + y_translation[t] - 1;
            }

            if (z_translation[t] >= 0) {//index increase
                zs0 = 0;
                ze0 = L - z_translation[t] - 1;
                zs1 = z_translation[t];
                ze1 = L - 1;
                //data1(xs1:xe1,ys1:ye1,zs1:ze1,t) = data1(xs0:xe0,ys0:ye0,zs0:ze0,t);
                for (int k = ze0, target = ze1; k >= zs0; --k, --target) {
                    cv::Rect srcRect(ys0, xs0, ye0 - ys0 + 1, xe0 - xs0 + 1);
                    cv::Rect dstRect(ys1, xs1, ye1 - ys1 + 1, xe1 - xs1 + 1);
                    data1[t][target](dstRect) = data1[t][k](srcRect).clone();
                }//for(k)
            } else {//index decrease
                zs0 = -z_translation[t];
                ze0 = L - 1;
                zs1 = 0;
                ze1 = L + z_translation[t] - 1;
                //data1(xs1:xe1,ys1:ye1,zs1:ze1,t) = data1(xs0:xe0,ys0:ye0,zs0:ze0,t);
                for (int k = zs0, target = zs1; k <= ze0; ++k, ++target) {
                    cv::Rect srcRect(ys0, xs0, ye0 - ys0 + 1, xe0 - xs0 + 1);
                    cv::Rect dstRect(ys1, xs1, ye1 - ys1 + 1, xe1 - xs1 + 1);
                    data1[t][target](dstRect) = data1[t][k](srcRect).clone();
                }//for(k)
            }//else
        }//for(t)

//      data1 = data1(max(x_translation)+1:end+min(x_translation),max(y_translation)+1:end+min(y_translation),
//      max(z_translation)+1:end+min(z_translation),:);
//        for (int i = 0; i < T; ++i) {
//            std::cout<<"x_translation["<<i<<"]: "<<x_translation[i]<<std::endl;
//        }
        int x_max = *std::max_element(x_translation, x_translation+T);
        int x_min = *std::min_element(x_translation, x_translation+T);
        int y_max = *std::max_element(y_translation, y_translation+T);
        int y_min = *std::min_element(y_translation, y_translation+T);
        int z_max = *std::max_element(z_translation, z_translation+T);
        int z_min = *std::min_element(z_translation, z_translation+T);
        if (x_min > 0){
            x_min = 0;
        }
        if (y_min > 0){
            y_min = 0;
        }
        if (z_min > 0){
            z_min = 0;
        }
//        std::cout<<"x_max:"<<x_max<<std::endl;
//        std::cout<<"x_min:"<<x_min<<std::endl;
//        std::cout<<"y_max:"<<y_max<<std::endl;
//        std::cout<<"y_min:"<<y_min<<std::endl;
//        std::cout<<"z_max:"<<z_max<<std::endl;
//        std::cout<<"z_min:"<<z_min<<std::endl;
        xs0 = x_max, xe0 = H - 1 + x_min;
        ys0 = y_max, ye0 = W - 1 + y_min;
        zs0 = z_max, ze0 = L - 1 + z_min;

        std::vector<std::vector<cv::Mat>> data_new(T);
        for (int t = 0; t < T; ++t) {
            for (int k = zs0; k<=ze0; ++k) {
//                std::cout<<"zs0: "<<zs0<<std::endl;
//                std::cout<<"ze0: "<<ze0<<std::endl;
//                if (data1[t][k].empty()) {
//                    std::cout << "empty" << std::endl;
//                    break;
//                }
                cv::Rect roi(ys0, xs0, ye0-ys0+1, xe0-xs0+1);
                data_new[t].emplace_back(data1[t][k](roi).clone());
            }
        }
//        for (int t = 0; t < T; ++t) {
//            for (int k = zs0; k <= ze0; ++k) {
////                std::cout<<data1[t][k].size()<<std::endl;
//                cv::Range rowRange(xs0, xe0);
//                cv::Range colRange(ys0, ye0);
//                data_new[t].emplace_back(data1[t][k](rowRange, colRange).clone());
//            }//for(k)
//        }//for(t)
        release3dMatrix(ref, H, W);
        release3dMatrix(moving,H,W);
        release3dMatrix(a_add,H_ext,W_ext);
        release3dMatrix(b_add,H_ext,W_ext);
        release3dMatrix(b_flip,H,W);
        release3dMatrix(matrix,H_ext,W_ext);
        delete[] array_1d;
//        delete[] matrix_1d;
        delete[] x_translation;
        delete[] y_translation;
        delete[] z_translation;

//        std::cout<<"height of image:"<< data_new[0][0].rows << std::endl;
//        std::cout<<"width of image:"<< data_new[0][0].cols  << std::endl;
//        std::cout<<"length of image:"<< data_new[0].size() << std::endl;
        return data_new;
    }//regCrossCorrelation()



}// namespace


