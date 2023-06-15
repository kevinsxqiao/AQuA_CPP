#include <iostream>
#include <boost/math/distributions/normal.hpp>

int main() {
    // 创建一个标准正态分布（均值为0，标准差为1）
    boost::math::normal_distribution<> standardNormal(0.0, 1.0);

    // 测试点
    double x = 1.0;

    // 计算概率密度函数(PDF)在x处的值
    double pdfValue = boost::math::pdf(standardNormal, x);

    // 计算累积分布函数(CDF)在x处的值
    double cdfValue = boost::math::cdf(standardNormal, x);

    // 输出结果
    std::cout << "PDF value at " << x << " is " << pdfValue << std::endl;
    std::cout << "CDF value at " << x << " is " << cdfValue << std::endl;

    return 0;
}
