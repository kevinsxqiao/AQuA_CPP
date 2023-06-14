#include <opencv2/opencv.hpp>
#include <iostream>

int main() {
    cv::Mat mat = (cv::Mat_<double>(3, 3) << 1, 2, 3,
            4, 5, 6,
            7, 8, 9);

    double minVal, maxVal;
    cv::Point minLoc, maxLoc;

    cv::minMaxLoc(mat, &minVal, &maxVal, &minLoc, &maxLoc);

    std::cout << "Min value: " << minVal << " at (" << minLoc.x << ", " << minLoc.y << ")\n";
    std::cout << "Max value: " << maxVal << " at (" << maxLoc.x << ", " << maxLoc.y << ")\n";

    return 0;
}
