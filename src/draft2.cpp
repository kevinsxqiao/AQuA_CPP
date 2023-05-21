#include <iostream>
#include <fftw3.h>

int main() {
    const int N = 8;

    double* in = (double*) fftw_malloc(sizeof(double) * N);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    fftw_plan plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

    for (int i = 0; i < N; ++i) {
        in[i] = i + 1;
    }

    fftw_execute(plan);

    std::cout << "result:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << out[i][0] << " + " << out[i][1] << "i" << std::endl;
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return 0;
}

//int main(){
//    return 1;
//}
