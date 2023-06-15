#include <opencv2/opencv.hpp>
#include <boost/math/distributions/normal.hpp>
#include <vector>
#include <cmath>

std::vector<cv::Mat> truncated_kept_var(const std::vector<cv::Mat>& quantiles) {
    int L = quantiles.size();
    std::vector<cv::Mat> pars(L);
    boost::math::normal_distribution<float> normal_dist; // standard normal distribution

    for (int k = 0; k < L; ++k) {
        const cv::Mat& q = quantiles[k];
        cv::Mat& p = pars[k];
        p = cv::Mat(q.rows, q.cols, CV_32F);

        for (int i = 0; i < q.rows; ++i) {
            for (int j = 0; j < q.cols; ++j) {
                float quantile = q.at<float>(i, j);

                if (quantile == 0.0f) {
                    p.at<float>(i, j) = 2.0f;
                } else {
                    float a = boost::math::quantile(normal_dist, quantile); // Equivalent to norminv in MATLAB
                    float phi_a = boost::math::pdf(normal_dist, a); // Equivalent to normpdf in MATLAB
                    float mu = a * quantile + phi_a;
                    float second_order = a * a * quantile + 1 - quantile + a * phi_a;
                    p.at<float>(i, j) = (second_order - mu * mu) * 2;
                }
            }
        }
    }

    return pars;
}

int main() {
    // Example usage:

    // Parameters
    int L = 2; // number of matrices
    int H = 3; // height
    int W = 3; // width

    // Creating the quantiles vector
    std::vector<cv::Mat> quantiles(L);
    for (int k = 0; k < L; ++k) {
        quantiles[k] = cv::Mat(H, W, CV_32F);
        // fill quantiles[k] with some values
    }

    // Compute pars
    std::vector<cv::Mat> pars = truncated_kept_var(quantiles);

    // Do something with pars
}
