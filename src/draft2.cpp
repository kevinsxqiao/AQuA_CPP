//gaussian filter
int dist = ceil(2 * opts.smoXY);
vector<cv::Mat> filter0(dist*2+1);
vector<cv::Mat> filter(dist*2+1);
//#pragma omp parallel for
for (int k = 0; k < dist * 2 + 1; ++k) {
filter0[k] = cv::Mat::zeros(dist*2+1, dist*2+1, CV_32F);
}
filter0[dist].at<float>(dist, dist) = 1;
int ksize = 2 * ceil(2 * opts.smoXY) + 1;
//#pragma omp parallel for
for (int k = 0; k < dist * 2 + 1; ++k) {
cv::GaussianBlur(filter0[k],filter0[k], cv::Size(ksize, ksize), opts.smoXY, opts.smoXY);
cv::pow(filter0[k], 2, filter[k]);
}

% Gaussian filter
dist = ceil(2*smo);
filter0 = zeros(dist*2+1,dist*2+1,dist*2+1);
filter0(dist+1,dist+1,dist+1) = 1;
filter0 = imgaussfilt(filter0,smo);
filter = filter0.^2;