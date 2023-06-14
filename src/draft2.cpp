#include <opencv2/opencv.hpp>
#include <vector>
#include <algorithm>
#include <numeric>

void baselineLinearEstimate(std::vector<std::vector<cv::Mat>>& data, int cut, int movAvgWin) {
    int T = data.size(); // 时间维度
    int L = data[0].size(); // 切片维度
    int H = data[0][0].rows; // 高度
    int W = data[0][0].cols; // 宽度

    // 计算移动平均值
    std::vector<std::vector<cv::Mat>> datMA(T, std::vector<cv::Mat>(L));
    for (int t = 0; t < T; t++) {
        for (int k = 0; k < L; k++) {
            cv::blur(data[t][k], datMA[t][k], cv::Size(movAvgWin, movAvgWin));
        }
    }

    int step = std::round(0.5 * cut);
    double maxVal;
    cv::minMaxLoc(cv::Mat(data[0][0]), nullptr, &maxVal);

    int nSegments = std::max(1, static_cast<int>(std::ceil(static_cast<double>(T) / step)) - 1);
    std::vector<cv::Mat> minPosition(H * W * L, cv::Mat::zeros(1, nSegments, CV_32S));

    for (int k = 0; k < nSegments; k++) {
        int t0 = 1 + k * step;
        int t1 = std::min(T, t0 + cut);

        for (int t = t0; t < t1; t++) {
            for (int l = 0; l < L; l++) {
                for (int i = 0; i < H; i++) {
                    for (int j = 0; j < W; j++) {
                        float curValue = datMA[t][l].at<float>(i, j);
                        if (curValue < datMA[t0][l].at<float>(i, j)) {
                            minPosition[i * W * L + j * L + l].at<int>(k) = t;
                        }
                    }
                }
            }
        }
    }

    std::vector<std::vector<cv::Mat>> F0(T, std::vector<cv::Mat>(L, cv::Mat::zeros(H, W, CV_32F)));

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            for (int l = 0; l < L; l++) {
                std::vector<int> curP;
                for (int k = 0; k < nSegments; k++) {
                    int pos = minPosition[i * W * L + j * L + l].at<int>(k);
                    if (pos != 0) {
                        curP.push_back(pos);
                    }
                }

                if (!curP.empty()) {
                    std::sort(curP.begin(), curP.end());
                    curP.erase(std::unique(curP.begin(), curP.end()), curP.end());

                    float firstVal = datMA[curP.front()][l].at<float>(i, j);
                    float lastVal = datMA[curP.back()][l].at<float>(i, j);

                    F0[0][l].at<float>(i, j) = firstVal;
                    F0[T-1][l].at<float>(i, j) = lastVal;

                    for (int m = 1; m < curP.size(); m++) {
                        int p1 = curP[m - 1];
                        int p2 = curP[m];
                        float v1 = datMA[p1][l].at<float>(i, j);
                        float v2 = datMA[p2][l].at<float>(i, j);

                        for (int t = p1; t < p2; t++) {
                            float interpVal = v1 + (v2 - v1) * (t - p1) / (p2 - p1);
                            F0[t][l].at<float>(i, j) = interpVal;
                        }
                    }
                } else {
                    for (int t = 0; t < T; t++) {
                        F0[t][l].at<float>(i, j) = maxVal;
                    }
                }
            }
        }
    }

    data = F0; // 最后把结果赋值回 data 变量
}
