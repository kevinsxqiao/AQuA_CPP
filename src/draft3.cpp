#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <mat.h>
#include <vector>
#include <boost/math/distributions/normal.hpp>
#include "data/data.h"
#include <omp.h>




int main(){
    vector<double> L_left = AQuA::load_vector("C:\\Users\\Kevin Qiao\\Desktop\\AQuA_data\\test\\L_left.bin");
    vector<double> L_right = AQuA::load_vector("C:\\Users\\Kevin Qiao\\Desktop\\AQuA_data\\test\\L_right.bin");
    vector<double> mus_Left = AQuA::load_vector("C:\\Users\\Kevin Qiao\\Desktop\\AQuA_data\\test\\mus_Left.bin");
    vector<double> mus_Right = AQuA::load_vector("C:\\Users\\Kevin Qiao\\Desktop\\AQuA_data\\test\\mus_Right.bin");
    double degreeOfFreedom_val = 710.0;
    vector<double>z_Left(L_left.size());
    vector<double>z_Right(L_right.size());
    //                !!!!!!!!!!!!!!!!LAST CORRECT!!!!!!!!!!!!!!!!!!!!!
    boost::math::normal_distribution<double> norm;
    for (int i = 0; i < L_left.size(); ++i) {
//        std::cout << "Iteration " << i << ": degreeOfFreedom_val=" << round(degreeOfFreedom_val) << ", mus_Left[i]=" << mus_Left[i] << std::endl;
        boost::math::non_central_t_distribution<double> nct(round(degreeOfFreedom_val), mus_Left[i]);
        double upperCDF = 1.0 - boost::math::cdf(nct, L_left[i]);
            if (upperCDF > 0.999999999999) {
                upperCDF = 0.999999999999;
            } else if (upperCDF < 0.000000000001) {
                upperCDF = 0.000000000001;
            }
//        cout<<"upperCDF:"<<upperCDF<<endl;
        z_Left[i] = -boost::math::quantile(norm, upperCDF);
    }
    for (int i = 0; i < L_right.size(); ++i) {
        boost::math::non_central_t_distribution<double> nct(round(degreeOfFreedom_val), mus_Right[i]);
        double upperCDF = 1.0 - boost::math::cdf(nct, L_right[i]);
        if (upperCDF > 0.999999999999) {
            upperCDF = 0.999999999999;
        } else if (upperCDF < 0.000000000001) {
            upperCDF = 0.000000000001;
        }
        z_Right[i] = -boost::math::quantile(norm, upperCDF);
    }
    return 0;
}