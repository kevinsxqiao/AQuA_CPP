#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/normal.hpp>

int main() {
    // 示例数据
    std::vector<double> L_left = {...}; // 你的输入数据
    std::vector<double> degreeOfFreedom = {...}; // 假设每个元素都有一个对应的自由度
    std::vector<double> mus_Left = {...}; // 你的输入数据

    int n = L_left.size(); // 假设所有数组都有相同的大小
    std::vector<double> z_Left(n);

    boost::math::normal_distribution<> norm;

    for (int i = 0; i < n; ++i) {
        // 计算非中心t分布的CDF的上尾值
        boost::math::non_central_t_distribution<> nct(round(degreeOfFreedom[i]), mus_Left[i]);
        double upperCDF = 1.0 - boost::math::cdf(nct, L_left[i]);

        // 计算正态分布的逆CDF
        z_Left[i] = boost::math::quantile(norm, upperCDF);
    }

    // 输出z_Left的值
    for (double z : z_Left) {
        std::cout << z << std::endl;
    }

    return 0;
}
