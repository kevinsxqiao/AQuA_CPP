//
// Created by Kevin Qiao on 2/17/2023.
//

#include "regCrossCorrelation.h"

namespace AQuA{


    DATA_TYPE medianFunc(DATA_TYPE* array, int size){
        DATA_TYPE median = 0;
        std::sort(array, array+size);
        if(size % 2 == 0){
            median = ( array[size/2] + array[size/2 - 1] ) / 2;
        }
        else{
            median = array[size/ 2];
        }
//        std::cout<< median;
        return median;
    }// medianFunc()


    float*** dft(float*** a_add, float*** b_add){// **need delete outside**  return to calCC(), release matrix in regCrossCorrelation()
//             float(&c)[H_ext][W_ext][L_ext]){
//        clock_t start,end;
//        start = clock();

// allocate memory  --use recommend fftwf_malloc which allocating memory on the heap,
// fixed-size array is allocated on the stack, using large fixed-size arrays may cause stack overflow issues or exceed the available stack size.
        float* input_a = static_cast<float *>(fftwf_malloc(H_ext * W_ext * L_ext * sizeof(float)));
        float* input_b = static_cast<float *>(fftwf_malloc(H_ext * W_ext * L_ext * sizeof(float)));
        fftwf_complex* output_a = static_cast<fftwf_complex *>(fftwf_malloc(H_ext * W_ext * (L_ext/2 + 1) * sizeof(fftwf_complex)));
        fftwf_complex* output_b = static_cast<fftwf_complex *>(fftwf_malloc(H_ext * W_ext * (L_ext/2 + 1) * sizeof(fftwf_complex)));
        fftwf_complex* output_dotProduct = static_cast<fftwf_complex *>(fftwf_malloc(H_ext * W_ext * (L_ext/2 + 1) * sizeof(fftwf_complex)));
        float* restored = static_cast<float *>(fftwf_malloc(H_ext * W_ext * L_ext * sizeof(float)));


        // initialize the input matrix (with sub2ind)
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                for (int k = 0; k < L_ext; ++k) {
                    input_a[i * W_ext * L_ext + j * L_ext + k] = a_add[i][j][k];
                    input_b[i * W_ext * L_ext + j * L_ext + k] = b_add[i][j][j];
//                    std::cout<<input_a[i * W_ext * L_ext + j * L_ext + k]<<" ";
                }
            }
        }
//        std::cout<<std::endl;

//    for (int k = 0; k < L_ext; ++k) {
//        for (int i = 0; i < H_ext; ++i) {
//            for (int j = 0; j < W_ext; ++j) {
//                std::cout<<input_a[i * W_ext * L_ext + j * L_ext + k]<<" ";
//            }
//            std::cout<<std::endl;
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<<std::endl;

        // create fftw plan
        fftwf_plan fft_plan_a = fftwf_plan_dft_r2c_3d(H_ext, W_ext, L_ext, input_a, output_a, FFTW_ESTIMATE);
        fftwf_plan fft_plan_b = fftwf_plan_dft_r2c_3d(H_ext, W_ext, L_ext, input_b, output_b, FFTW_ESTIMATE);
        fftwf_plan ifft_plan = fftwf_plan_dft_c2r_3d(H_ext, W_ext, L_ext, output_dotProduct, restored,  FFTW_ESTIMATE);


        fftwf_execute(fft_plan_a);
        fftwf_execute(fft_plan_b);


//print the matrix after fft_plan_a
//    for (int k = 0; k < L_ext ; ++k) {
//        for (int j = 0; j < W_ext; ++j) {
//            for (int i = 0; i < H_ext; ++i) {
//                int index = i * W_ext * (L_ext/2 + 1) + j * (L_ext/2 + 1) + k;
//                std::cout<< output_a[index][0]<< "+" << output_a[index][1]<< "i"<< "      ";
//            }
//            std::cout<< std::endl;
//        }
//        std::cout<< std::endl;
//    }


        // dot product
        for (int index = 0; index<H_ext * W_ext * (L_ext/2 + 1); ++index) {
            output_dotProduct[index][0] = output_a[index][0] * output_b[index][0] - output_a[index][1] * output_b[index][1];
            output_dotProduct[index][1] = output_a[index][0] * output_b[index][1] + output_a[index][1] * output_b[index][0];
        }


        fftwf_execute(ifft_plan);

        float*** c = create3dMatrix_ext_float();
        // normalize the IFFT output, since ifftn() will automatically normalize in matlab
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                for (int k = 0; k < L_ext; ++k) {
                    restored[i * W_ext * L_ext + j * L_ext + k] /= H_ext*W_ext*L_ext;
                    c[i][j][k] = restored[i * W_ext * L_ext + j * L_ext + k];
                }
            }
        }

        // free
        fftwf_destroy_plan(fft_plan_a);
        fftwf_destroy_plan(fft_plan_b);
        fftwf_destroy_plan(ifft_plan);
        fftwf_free(input_a);
        fftwf_free(input_b);
        fftwf_free(output_a);
        fftwf_free(output_b);
        fftwf_free(restored);
        fftwf_free(output_dotProduct);

//        end=clock();
//        std::cout<< "time:"<< double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
        return c;
    }


// c = ifftn(fftn(a_add).*fftn(b_add));  //for double type
    double*** dft(double*** a_add, double*** b_add){// **need delete outside**  return to calCC(), release matrix in regCrossCorrelation()
//        clock_t start,end;
//        start = clock();

        // allocate memory  --use recommend fftw_malloc which allocating memory on the heap,
        // fixed-size array is allocated on the stack, using large fixed-size arrays may cause stack overflow issues or exceed the available stack size.
        double* input_a = static_cast<double *>(fftw_malloc(H_ext * W_ext * L_ext * sizeof(double)));
        double* input_b = static_cast<double *>(fftw_malloc(H_ext * W_ext * L_ext * sizeof(double)));
        fftw_complex* output_a = static_cast<fftw_complex *>(fftw_malloc(H_ext * W_ext * (L_ext/2 + 1) * sizeof(fftw_complex)));
        fftw_complex* output_b = static_cast<fftw_complex *>(fftw_malloc(H_ext * W_ext * (L_ext/2 + 1) * sizeof(fftw_complex)));
        fftw_complex* output_dotProduct = static_cast<fftw_complex *>(fftw_malloc(H_ext * W_ext * (L_ext/2 + 1) * sizeof(fftw_complex)));
        double* restored = static_cast<double *>(fftw_malloc(H_ext * W_ext * L_ext * sizeof(double)));


        // initialize the input matrix (with sub2ind)
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                for (int k = 0; k < L_ext; ++k) {
                    input_a[i * W_ext * L_ext + j * L_ext + k] = a_add[i][j][k];
                    input_b[i * W_ext * L_ext + j * L_ext + k] = b_add[i][j][j];
//                    std::cout<<input_a[i * W_ext * L_ext + j * L_ext + k]<<" ";
                }
            }
        }
//        std::cout<<std::endl;

//    for (int k = 0; k < L_ext; ++k) {
//        for (int i = 0; i < H_ext; ++i) {
//            for (int j = 0; j < W_ext; ++j) {
//                std::cout<<input_a[i * W_ext * L_ext + j * L_ext + k]<<" ";
//            }
//            std::cout<<std::endl;
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<<std::endl;

        // create fftw plan
        fftw_plan fft_plan_a = fftw_plan_dft_r2c_3d(H_ext, W_ext, L_ext, input_a, output_a, FFTW_ESTIMATE);
        fftw_plan fft_plan_b = fftw_plan_dft_r2c_3d(H_ext, W_ext, L_ext, input_b, output_b, FFTW_ESTIMATE);
        fftw_plan ifft_plan = fftw_plan_dft_c2r_3d(H_ext, W_ext, L_ext, output_dotProduct, restored,  FFTW_ESTIMATE);


        fftw_execute(fft_plan_a);
        fftw_execute(fft_plan_b);


//print the matrix after fft_plan_a
//    for (int k = 0; k < L_ext ; ++k) {
//        for (int j = 0; j < W_ext; ++j) {
//            for (int i = 0; i < H_ext; ++i) {
//                int index = i * W_ext * (L_ext/2 + 1) + j * (L_ext/2 + 1) + k;
//                std::cout<< output_a[index][0]<< "+" << output_a[index][1]<< "i"<< "      ";
//            }
//            std::cout<< std::endl;
//        }
//        std::cout<< std::endl;
//    }


        // dot product
        for (int index = 0; index<H_ext * W_ext * (L_ext/2 + 1); ++index) {
            output_dotProduct[index][0] = output_a[index][0] * output_b[index][0] - output_a[index][1] * output_b[index][1];
            output_dotProduct[index][1] = output_a[index][0] * output_b[index][1] + output_a[index][1] * output_b[index][0];
        }


        fftw_execute(ifft_plan);

        double*** c = create3dMatrix_ext();
        // normalize the IFFT output, since ifftn() will automatically normalize in matlab
        for (int i = 0; i < H_ext; ++i) {
            for (int j = 0; j < W_ext; ++j) {
                for (int k = 0; k < L_ext; ++k) {
                    restored[i * W_ext * L_ext + j * L_ext + k] /= H_ext*W_ext*L_ext;
                    c[i][j][k] = restored[i * W_ext * L_ext + j * L_ext + k];
                }
            }
        }

        // free
        fftw_destroy_plan(fft_plan_a);
        fftw_destroy_plan(fft_plan_b);
        fftw_destroy_plan(ifft_plan);
        fftw_free(input_a);
        fftw_free(input_b);
        fftw_free(output_a);
        fftw_free(output_b);
        fftw_free(restored);
        fftw_free(output_dotProduct);

//        end=clock();
//        std::cout<< "time:"<< double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
        return c;
    }


    // 90-degree anticlockwise rotation of matrix
    DATA_TYPE*** rotate2d(DATA_TYPE*** ref){  // not release yet
        DATA_TYPE*** mat = create3dMatrix();
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                mat[H - j][i][0] = ref[i][j][0];
            }// for(j)
        }// for(i)
        return mat;
    }// rotate2dMatrix()



// flip(flip(flip(b,1),2),3);       **need delete outside**
//    void flip3d(DATA_TYPE (&ref)[H][W][L],
//                 DATA_TYPE (&mat)[H][W][L]){
    DATA_TYPE*** flip3d(DATA_TYPE*** ref, DATA_TYPE*** b_flip){


        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    b_flip[H - 1 - i][W - 1 -j ][L - 1 - k] = ref[i][j][k];
                }// for(k)
            }// for(j)
        }// for(i)

        return b_flip;
    }// flip_3d()

// input: a = moving; b = ref output: c[]
    DATA_TYPE*** calCC(DATA_TYPE*** a, DATA_TYPE*** b, DATA_TYPE*** a_add, DATA_TYPE*** b_add, DATA_TYPE*** b_flip){

//        DATA_TYPE b_flipped[H][W][L];
//        DATA_TYPE a_add[H_ext][W_ext][L_ext];
//        DATA_TYPE b_add[H_ext][W_ext][L_ext];
        b_flip = flip3d(b,b_flip); // flip(flip(flip(b,1),2),3); **need delete outside**
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



    DATA_TYPE**** regCrossCorrelation(DATA_TYPE**** data1, DATA_TYPE**** data2){
        /*
         * create variables and allocate memory
         */

//        int x_translation[T], y_translation[T], z_translation[T];
//        DATA_TYPE ref[H][W][L];
//        DATA_TYPE temp_col[H * W * L];// convert moving and ref to one-dimension for sorting
//        DATA_TYPE moving[H][W][L];
//        DATA_TYPE matrix[H_ext][W_ext][L_ext];// store the output of calCC()
//        DATA_TYPE matrix_col[H_ext * W_ext * L_ext];
        DATA_TYPE sum, median;
        int frame_start = 0;
        int frame_end = 2;
        int wShift, hShift, lShift;
        int id,r1;
        int xs0, xe0, xs1, xe1, ys0, ye0, ys1, ye1, zs0, ze0, zs1, ze1;
        int* x_translation = new int [T];
        int* y_translation = new int [T];
        int* z_translation = new int [T];
        DATA_TYPE*** ref = create3dMatrix();
        DATA_TYPE* temp_col = new DATA_TYPE [H*W*L];
        DATA_TYPE*** moving = create3dMatrix();
        DATA_TYPE*** a_add = create3dMatrix_ext();
        DATA_TYPE*** b_add = create3dMatrix_ext();
        DATA_TYPE*** b_flip = create3dMatrix();
        DATA_TYPE*** matrix = create3dMatrix();
        DATA_TYPE* matrix_col = new DATA_TYPE [H_ext*W_ext*L_ext];
        DATA_TYPE ****dat1;// released in main.cpp: releaseData(datOrg1);



        // compute the mean value from frame_start to frame_end

        for (int i=0, m=0; i<H;++i) {
            for (int j=0; j<W;++j) {
                for(int k=0; k<L;++k){
                    sum = 0;
                    for (int t=frame_start;t<frame_end;++t) {
                        sum += data1[i][j][k][t];
                    }// for(t)
                    ref[i][j][k] = sum / DATA_TYPE(frame_end - frame_start);
                    temp_col[m++] = ref[i][j][k]; //convert ref to one dimension
                }// for(k)
            }// for(j)
        }// for(i)

        // compute the median value of temp_col --- medianFunc(<array>, <size>)
        median = medianFunc(temp_col, H * W * L);


        // align dark
        for(int i=0;i<H;++i){
            for(int j=0;j<W;++j){
                for(int k=0;k<L;++k){
                    ref[i][j][k] -= median;  // ref = ref - median(ref(:));
                }// for(k)
            }// for(j)
        }// for(i)


        // concise cross correlation --loop from frame <2> to <T>
        for(int t=0; t<T; ++t){
            for(int i=0, m=0; i<H; ++i){
                for(int j=0; j<W; ++j){
                    for(int k=0; k<L; ++k){
                        moving[i][j][k] = data1[i][j][k][t]; //initialize moving
                        temp_col[m++] = moving[i][j][k]; //convert moving to one dimension
                    }// for(k)
                } // for(j)
            }// for(i)
            median = medianFunc(temp_col, H * W * L);
            for(int i=0;i<H;++i){
                for(int j=0;j<W;++j){
                    for(int k=0;k<L;++k){
                        moving[i][j][k] -= median;  // moving = moving - median(moving(:));
                        matrix = calCC(moving, ref, a_add, b_add, b_flip); // matrix = calCC(moving,ref);
                    }// for(k)
                } //for(j)
            }// for(i)

            for(int i=0, m=0;i<H_ext;++i){
                for(int j=0;j<W_ext;++j){
                    for(int k=0;k<L_ext;++k){
                        matrix_col[m] = matrix[i][j][k];  // convert matrix to one-dimension to find the position of the maximum element
                    }// for(k)
                } //for(j)
            }// for(i)

            id = std::max_element(matrix_col, matrix_col+(H_ext*W_ext*L_ext)) - matrix_col;  //[~,id] = max(matrix(:));
            // [hShift,wShift,lShift] = ind2sub(size(matrix),id);
            hShift = id/(W_ext*L_ext); // x as page, in which page
            r1 = id%(W_ext*L_ext); // in the page hShift(2d matrix), linear position
            wShift = r1 / L_ext; // in the page hShift(2d matrix), which row;
            lShift = r1 % L_ext; // in the page hShift(2d matrix), which column;
            x_translation[t] = H - hShift;
            y_translation[t] = W - wShift;
            z_translation[t] = L - lShift;
        }// for(t)


        for(int t=0;t<T;++t){
            if (x_translation[t]>=0){
                xs0 = 1;
                xe0 = H - x_translation[t];
                xs1 = 1 + x_translation[t];
                xe1 = H;
            }
            else{
                xs0 = 1 - x_translation[t];
                xe0 = H;
                xs1 = 1;
                xe1 = H + x_translation[t];
            }

            if (y_translation[t]>=0){
                ys0 = 1;
                ye0 = W - y_translation[t];
                ys1 = 1 + y_translation[t];
                ye1 = W;
            }
            else{
                ys0 = 1 - y_translation[t];
                ye0 = W;
                ys1 = 1;
                ye1 = W + y_translation[t];
            }

            if (z_translation[t]>=0){
                zs0 = 1;
                ze0 = L - z_translation[t];
                zs1 = 1 + z_translation[t];
                ze1 = L;
            }
            else{
                zs0 = 1 - z_translation[t];
                ze0 = L;
                zs1 = 1;
                ze1 = L + z_translation[t];
            }

            //data1(xs1:xe1,ys1:ye1,zs1:ze1,t) = data1(xs0:xe0,ys0:ye0,zs0:ze0,t);
            for (int i = xs0, x_target = xs1; i <= xe0; ++i, ++x_target) {
                for (int j = ys0, y_target = ys1; j <= ye0; ++j, ++y_target) {
                    for (int k = zs0, z_target = zs1; k <= zs1; ++k, ++z_target) {
                        data1[x_target][y_target][z_target][t] = data1[i][j][k][t];
                    } //for(k)
                } //for(j)
            } //for(i)
//            if(~isempty(data2))
//                data2(xs1:xe1,ys1:ye1,zs1:ze1,t) = data2(xs0:xe0,ys0:ye0,zs0:ze0,t);
//            end
        } //for(t)
        int x_max = std::max_element(x_translation, x_translation+T) - x_translation;
        int x_min = std::min_element(x_translation, x_translation+T) - x_translation;
        int y_max = std::max_element(y_translation, y_translation+T) - y_translation;
        int y_min = std::min_element(y_translation, y_translation+T) - y_translation;
        int z_max = std::max_element(z_translation, z_translation+T) - z_translation;
        int z_min = std::min_element(z_translation, z_translation+T) - z_translation;


        int x = x_max - x_min + 1, y = y_max - y_min + 1, z = z_max - z_min+ 1; // x,y,z indicates the size of the new matrix after registration
        dat1 = new DATA_TYPE ***[x];
        for (int i = 0; i < x; ++i) {
            dat1[i] = new DATA_TYPE**[y];

        }
        for (int i = 0; i < x; ++i) {
            for (int j = 0; j < y; ++j) {
                dat1[i][j] = new DATA_TYPE*[z];
            }
        }
        for (int i = 0; i < x; ++i) {
            for (int j = 0; j < y; ++j) {
                for (int k = 0; k < z; ++k) {
                    dat1[i][j][k] = new DATA_TYPE[T];
                }
            }
        }
//        data1 = data1(max(x_translation)+1:end+min(x_translation),max(y_translation)+1:end+min(y_translation),max(z_translation)+1:end+min(z_translation),:);
        for (int i = x_min, x_target=0; i < x_max+1; ++i, ++x_target) {
            for (int j =y_min, y_target=0; j < y_max+1; ++j, ++y_target) {
                for (int k = z_min, z_target=0; k < z_max+1; ++k, ++z_target) {
                    for (int t = 0; t < T; ++t) {
                        dat1[x_target][y_target][z_target][t] = data1[i][j][k][t];
                    }
                }
            }
        }
//        if(~isempty(data2))
//            data2 = data2(max(x_translation)+1:end+min(x_translation),max(y_translation)+1:end+min(y_translation),max(z_translation)+1:end+min(z_translation),:);
//        end

        releaseData(x_translation);
        releaseData(y_translation);
        releaseData(z_translation);
        releaseData(ref);
        releaseData(temp_col);
        releaseData(moving);
        releaseData(a_add);
        releaseData(b_add);
        releaseData(b_flip);
        releaseData(matrix);
        releaseData(matrix_col);

        return dat1;
    }// regCrossCorrelation()



}// namespace


