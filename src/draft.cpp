#include <iostream>
#include <vector>
#include <cmath> // For std::round()
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/normal.hpp>

int main() {
    // Sample data
    std::vector<float> L_left = {1.5, 2.0, -2.5};
    double degreeOfFreedomValue = 10.0;
    std::vector<float> mus_Left = {0.5, 1.0, -1.5};

    // Ensure all vectors have the same size, otherwise handle error
    int n = L_left.size();
    std::vector<float> z_Left(n);

    for (int i = 0; i < n; ++i) {
        // Calculate the upper tail CDF value of the non-central t-distribution
        boost::math::non_central_t_distribution<float> nct(std::round(degreeOfFreedomValue), mus_Left[i]);
        double upperCDF = 1.0 - boost::math::cdf(nct, L_left[i]);

        // Calculate the inverse CDF of the normal distribution
        boost::math::normal_distribution<float> norm;
        z_Left[i] = -boost::math::quantile(norm, upperCDF);
    }

    // Output z_Left values
    for (double z : z_Left) {
        std::cout << z << std::endl;
    }

    return 0;
}
