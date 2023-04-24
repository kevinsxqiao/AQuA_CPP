#include <opencv2/opencv.hpp>
#include <iostream>

int main() {
    std::vector<cv::Mat> volumes;
    cv::imreadmulti("C:/Users/Kevin Qiao/Desktop/3D_data/3D_dataFrame 1.tif", volumes, cv::IMREAD_ANYDEPTH);

    // 检查图像是否成功加载
    if (volumes.empty()) {
        std::cerr << "Failed to read volume file!" << std::endl;
        return 1;
    }

    // 输出体积的大小和数据类型
    cv::Mat volume = volumes[0];
    std::cout << "Volume size: " << volume.size() << std::endl;
    std::cout << "Volume type: " << volume.type() << std::endl;

    return 0;
}
