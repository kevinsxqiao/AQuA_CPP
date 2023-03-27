#include <iostream>
#include <cmath>
#include <fftw3.h>
#include<time.h>

#define N0 2
#define N1 2
#define N2 2
#define N0_ext (2*N0-1)
#define N1_ext (2*N1-1)
#define N2_ext (2*N2-1)

void dft(double a_add[N0_ext][N1_ext][N2_ext], double b_add[N0_ext][N1_ext][N2_ext], double(&c)[N0_ext][N1_ext][N2_ext]){
    clock_t start,end;
    start = clock();

    // allocate memory  --use recommend fftw_malloc which allocating memory on the heap,
    // fixed-size array is allocated on the stack, using large fixed-size arrays may cause stack overflow issues or exceed the available stack size.
    double* input_a = static_cast<double *>(fftw_malloc(N0_ext * N1_ext * N2_ext * sizeof(double)));
    double* input_b = static_cast<double *>(fftw_malloc(N0_ext * N1_ext * N2_ext * sizeof(double)));
    fftw_complex* output_a = static_cast<fftw_complex *>(fftw_malloc(N0_ext * N1_ext * (N2_ext/2 + 1) * sizeof(fftw_complex)));
    fftw_complex* output_b = static_cast<fftw_complex *>(fftw_malloc(N0_ext * N1_ext * (N2_ext/2 + 1) * sizeof(fftw_complex)));
    fftw_complex* output_dotProduct = static_cast<fftw_complex *>(fftw_malloc(N0_ext * N1_ext * (N2_ext/2 + 1) * sizeof(fftw_complex)));
    double* restored = static_cast<double *>(fftw_malloc(N0_ext * N1_ext * N2_ext * sizeof(double)));


    // initialize the input matrix (with sub2ind)
    for (int i = 0; i < N0_ext; ++i) {
        for (int j = 0; j < N1_ext; ++j) {
            for (int k = 0; k < N2_ext; ++k) {
                input_a[i * N1_ext * N2_ext + j * N2_ext + k] = a_add[i][j][k];
                input_b[i * N1_ext * N2_ext + j * N2_ext + k] = b_add[i][j][j];
                std::cout<<input_a[i * N1_ext * N2_ext + j * N2_ext + k]<<" ";
            }
        }
    }
    std::cout<<std::endl;

//    for (int k = 0; k < N2_ext; ++k) {
//        for (int i = 0; i < N0_ext; ++i) {
//            for (int j = 0; j < N1_ext; ++j) {
//                std::cout<<input_a[i * N1_ext * N2_ext + j * N2_ext + k]<<" ";
//            }
//            std::cout<<std::endl;
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<<std::endl;

    // create fftw plan
    fftw_plan fft_plan_a = fftw_plan_dft_r2c_3d(N0_ext, N1_ext, N2_ext, input_a, output_a, FFTW_ESTIMATE);
    fftw_plan fft_plan_b = fftw_plan_dft_r2c_3d(N0_ext, N1_ext, N2_ext, input_b, output_b, FFTW_ESTIMATE);
    fftw_plan ifft_plan = fftw_plan_dft_c2r_3d(N0_ext, N1_ext, N2_ext, output_dotProduct, restored,  FFTW_ESTIMATE);


    fftw_execute(fft_plan_a);
    fftw_execute(fft_plan_b);


//print the matrix after fft_plan_a
//    for (int k = 0; k < N2_ext ; ++k) {
//        for (int j = 0; j < N1_ext; ++j) {
//            for (int i = 0; i < N0_ext; ++i) {
//                int index = i * N1_ext * (N2_ext/2 + 1) + j * (N2_ext/2 + 1) + k;
//                std::cout<< output_a[index][0]<< "+" << output_a[index][1]<< "i"<< "      ";
//            }
//            std::cout<< std::endl;
//        }
//        std::cout<< std::endl;
//    }


    // dot product
    for (int index = 0; index<N0_ext * N1_ext * (N2_ext/2 + 1); ++index) {
                output_dotProduct[index][0] = output_a[index][0] * output_b[index][0] - output_a[index][1] * output_b[index][1];
                output_dotProduct[index][1] = output_a[index][0] * output_b[index][1] + output_a[index][1] * output_b[index][0];
    }


    fftw_execute(ifft_plan);


    // normalize the IFFT output, since ifftn() will automatically normalize in matlab
    for (int i = 0; i < N0_ext; ++i) {
        for (int j = 0; j < N1_ext; ++j) {
            for (int k = 0; k < N2_ext; ++k) {
                restored[i * N1_ext * N2_ext + j * N2_ext + k] /= N0_ext*N1_ext*N2_ext;
                c[i][j][k] = restored[i * N1_ext * N2_ext + j * N2_ext + k];
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

    end=clock();
    std::cout<< "time:"<< double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
}

int main(){
    double a[N0_ext][N1_ext][N2_ext]={1,2,3,4,55,3,3,6,2,3,4,4,6,4,4,43,3,24,54,54,52,45,245,2345,45,5,45};
    double b[N0_ext][N1_ext][N2_ext]={32,3,3,3,7,5,3,2,6,7,8,5,3,6,7,5,32,2,4,6,7,3,4,2,5,6,5.2345345};
    double c[N0_ext][N1_ext][N2_ext]= {0};
    dft(a,b,c);
    for (int k = 0; k < N2_ext; ++k) {
        for (int i = 0; i < N0_ext; ++i) {
            for (int j = 0; j < N1_ext; ++j) {
                std::cout<< c[i][j][k]<< " ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
    return 0;
}
