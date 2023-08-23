//
// Created by Kevin Qiao on 3/28/2023.
//

#include "baselineRemoveAndNoiseEstimation.h"


namespace AQuA{


    vector<cv::Mat> fit_F0_var(const vector<cv::Mat>& F0ProOrg, const vector<cv::Mat>& varMapOrg, int H, int W, int L, int T, int dist) {
        vector<cv::Mat> varMapOut((vector<cv::Mat>(L)));
        vector<cv::Mat> F0Pro ((vector<cv::Mat>(L)));
        vector<cv::Mat> varMap ((vector<cv::Mat>(L)));
//        vector<float> F0Pro;
//        vector<float> varMap;
//        vector<float> select1;
        vector<cv::Mat> select1 ((vector<cv::Mat>(L)));
        vector<cv::Mat> select2 ((vector<cv::Mat>(L)));
        
//#pragma omp parallel for
        for (int k = 0; k < L; ++k) {
            F0Pro[k] = cv::Mat(H-2*dist,W-2*dist,CV_32F);
            varMapOut[k] = cv::Mat(H,W,CV_32F);
            varMap[k] = cv::Mat(H-2*dist,W-2*dist,CV_32F);
        }

        //remove values close to 0
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            cv::Rect roi(dist,dist,H-dist,W-dist);
            F0Pro[k] = F0ProOrg[k](roi).clone();
            varMap[k] = varMapOrg[k](roi).clone();
            for (int i = dist; i < H-dist; ++i) {
                for (int j = dist; j < W-dist; ++j) {
                    if (varMapOrg[k].at<float>(i,j) < 1e-8) {
                        F0Pro[k].at<float>(i, j) = NAN;
                        varMap[k].at<float>(i, j) = NAN;
                    }
                }
            }
        }
//        for (int k = 0; k < L; ++k) {
//            for (int i = dist; i < H-dist; ++i) {
//                for (int j = dist; j < W-dist; ++j) {
//                    if (varMapOrg[k].at<float>(i,j) > 1e-8) {
//                        F0Pro.push_back(F0ProOrg[k].at<float>(i,j));
//                        varMap.push_back(varMapOrg[k].at<float>(i,j));
//                    }
//                }
//            }
//        }
//        cout<<"F0Pro: "<<endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<F0Pro[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
//        if isempty(F0Pro)
//           % if no valid variance
//        varMapOut = ones(size(F0ProOrg)) * 1e-8;
//        return;
//        end

        //downsample --- the number of original pixels is too large.
        float maxX = 0;
        float minX = 1;
        double maxVal;
        double minVal;
//        AQuA::writeDataToMatFile(F0Pro, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test_F0.mat");
//        cv::Point *maxLoc;
        for (int k = 0; k < L; ++k) {
            cv::minMaxLoc(cv::Mat(F0Pro[k]), &minVal, &maxVal);
            if (maxVal > maxX){
                maxX = static_cast<float>(maxVal);
            }
            if (minVal < minX){
                minX = static_cast<float>(minVal);
            }
        }
//        if (dist == 1){
//            cout<<"countV: "<<countV<<endl;
//            cout<<"min max value:"<< minX<< " "<< maxX<<endl;
//            cout<<"F0Pro: "<<endl;
//            for (int k = 0; k < 5; ++k) {
//                for (int i = 0; i < 20; ++i) {
//                    for (int j = 0; j < 20; ++j) {
//                        cout<<F0Pro[0].at<float>(i,j)<<" ";
//                    }
//                    cout<<endl;
//                }
//            }
//        }
        float delta = max(static_cast<float>(1e-5), (maxX - minX)/2000);
        cv::Mat x = cv::Mat::zeros(2000,1,CV_32F);
        cv::Mat y = cv::Mat::zeros(2000,1,CV_32F);
        cv::Mat valid = cv::Mat::zeros(2000,1,CV_32S);
        for (int ii = 0; ii < 2000; ++ii) {
            bool flag = false;
//#pragma omp parallel for collapse(3) shared(flag)
            for (int k = 0; k < L; ++k) {
                select1[k] = cv::Mat::zeros(H-2*dist,W-2*dist,CV_32S);
                for (int i = 0; i < H-2*dist; ++i) {
                    for (int j = 0; j < W-2*dist; ++j) {
                        if (ii == 0 ){
                            if ((F0Pro[k].at<float>(i,j) >= minX + ii*delta) && (F0Pro[k].at<float>(i,j) <= minX + (ii+1)*delta)){
                                select1[k].at<int>(i,j) = 1;
                                flag = true;
                            }
                        }else{
                            if ((F0Pro[k].at<float>(i,j) > minX + ii*delta) && (F0Pro[k].at<float>(i,j) < minX + (ii+1)*delta)){
                                select1[k].at<int>(i,j) = 1;
                                flag = true;
                            }
                        }
                    }//for(j)
                }//for(i)
            }//for(k)
            if(flag){
                float sumX = 0;
                float sumY = 0;
                int count = 0;
                for (int k = 0; k < L; ++k) {
                    for (int i = 0; i < H-2*dist; ++i) {
                        for (int j = 0; j < W-2*dist; ++j) {
                            if (select1[k].at<int>(i,j) == 1){
                                sumX += F0Pro[k].at<float>(i,j);
                                sumY += varMap[k].at<float>(i,j);
                                ++count;
                            }
                        }
                    }
                }
                x.at<float>(ii,0) = sumX / static_cast<float>(count);
                y.at<float>(ii,0) = sumY / static_cast<float>(count);
                valid.at<int>(ii,0) = 1;
            }//if(flag)
        }//for(ii)
        vector<float> x_new;
        vector<float> y_new;
        for (int i = 0; i < 2000; ++i) {
            if (valid.at<int>(i,0) == 1){
                x_new.push_back(x.at<float>(i,0));
                y_new.push_back(y.at<float>(i,0));
            }
        }

        if (x_new.size() == 1){
//#pragma omp parallel for
            for (int k = 0; k < L; ++k) {
                varMapOut[k] = y_new[0];
            }
            return varMapOut;
        }
        //graph construction --- the start point and end point, if pick the most extreme ones, too
        // unstable. Since source is denser, we could use more points.
        int source = ceil(y_new.size() * 0.05);
        int sink = max(static_cast<double>(1), floor(y_new.size() * 0.99));
        //first layer
        vector<double> dist1(y_new.size(), numeric_limits<double>::infinity());
        vector<int> preMap1(y_new.size(), 0);
////#pragma omp parallel for collapse(3)
        for (int j = 0; j < y_new.size(); ++j) {
            double minCost = numeric_limits<double>::infinity();
            int preNode = j;
            for (int i = 0; i < min(j - 1, source); ++i) {
                double a = (y_new[j] - y_new[i]) / (x_new[j] - x_new[i]);
                double b = y_new[i] - a * x_new[i];
//                % first term => first inclined segment
//                % second term => horizontal
                double cost = 0;
                double cost1 = 0;
                double cost2 = 0;
                for (int k = i; k < j; ++k) {
                    cost1 += abs(y_new[k] - a*x_new[k]+b);
                }
                for (int k = 0; k < i; ++k) {
                    cost2 += abs(y_new[k] - y_new[i]);
                }
                cost = cost1 + cost2;
                if (cost < minCost){
                    minCost = cost;
                    preNode = i;
                }
                dist1[j] = minCost;
                preMap1[j] = preNode;
            }//for(i)
        }//for(j)
        //second layer
        vector<double> dist2(y_new.size(), numeric_limits<double>::infinity());
        vector<int> preMap2(y_new.size(), 0);
////#pragma omp parallel for collapse(3)
        for (int j = 0; j < y_new.size(); ++j) {
            double minCost = numeric_limits<double>::infinity();
            int preNode = j;
            for (int i = 0; i < j-1; ++i) {
                double a = (y_new[j] - y_new[i]) / (x_new[j] - x_new[i]);
                double b = y_new[i] - a * x_new[i];
                double cost = 0;
                for (int k = i; k < j; ++k) {
                    cost += abs(y_new[k] - a*x_new[k]+b);
                }
                if (dist1[i] + cost < minCost){
                    minCost = dist1[i] + cost;
                    preNode = i;
                }
                dist2[j] = minCost;
                preMap2[j] = preNode;
            }//for(i)
        }//for(j)
        //third layer
        vector<double> dist3(y_new.size(), numeric_limits<double>::infinity());
        vector<int> preMap3(y_new.size(), 0);
////#pragma omp parallel for collapse(3)
        for (int j = sink-1; j < y_new.size(); ++j) {
            double minCost = numeric_limits<double>::infinity();
            int preNode = j;
            for (int i = 0; i < j-1; ++i) {
                double a = (y_new[j] - y_new[i]) / (x_new[j] - x_new[i]);
                double b = y_new[i] - a * x_new[i];
                double cost = 0;
                double cost1 = 0;
                double cost2 = 0;
                for (int k = i; k < j; ++k) {
                    cost1 += abs(y_new[k] - a*x_new[k]+b);
                }
                for (int k = j; k < y_new.size(); ++k) {
                    cost2 += abs(y_new[k] - y_new[j]);
                }
                cost = cost1 + cost2;
                if (dist2[i] + cost < minCost){
                    minCost = dist2[i] + cost;
                    preNode = i;
                }
                dist3[j] = minCost;
                preMap3[j] = preNode;
            }//for(i)
        }//for(j)
        //sink
        int node3 = min_element(dist3.begin(), dist3.end()) - dist3.begin();
        double x3 = x_new[node3];
        double y3 = y_new[node3];
        //end of 2nd segment
        int node2 = preMap3[node3];
        double x2 = x_new[node2];
        double y2 = y_new[node2];
        //end of 1st segment
        int node1 = preMap2[node2];
        double x1 = x_new[node1];
        double y1 = y_new[node1];
        //start of 1st segment
        int node0 = preMap1[node1];
        double x0 = x_new[node0];
        double y0 = y_new[node0];

        double a1 = (y0 - y1) / (x0 -x1);
        double b1 = y1 - a1*x1;

        double a2 = (y1 - y2) / (x1 -x2);
        double b2 = y2 - a2*x2;

        double a3 = (y3 - y2) / (x3 -x2);
        double b3 = y3 - a3*x3;
        //fitting
//#pragma omp parallel for
        for (int k = 0; k < L; ++k) {
            varMapOut[k] = varMapOrg[k];
        }
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            select2[k] = cv::Mat(H,W,CV_32S);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    if (F0ProOrg[k].at<float>(i,j) <=x0){
                        varMapOut[k].at<float>(i,j) = static_cast<float>(y0);
                    }
                    if ((F0ProOrg[k].at<float>(i,j) >= x0) && (F0ProOrg[k].at<float>(i,j) < x1)){
                        select2[k].at<int>(i,j) = 1;
                    } else{
                        select2[k].at<int>(i,j) = 0;
                    }
                }
            }
        }
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    if (select2[k].at<int>(i,j) ==1){
                        varMapOut[k].at<float>(i,j) = static_cast<float>(a1 * F0ProOrg[k].at<float>(i,j) + b1);
                    }
                    if ((F0ProOrg[k].at<float>(i,j) >= x1) && (F0ProOrg[k].at<float>(i,j) < x2)){
                        select2[k].at<int>(i,j) = 1;
                    } else{
                        select2[k].at<int>(i,j) = 0;
                    }
                }
            }
        }
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    if (select2[k].at<int>(i,j) ==1){
                        varMapOut[k].at<float>(i,j) = static_cast<float>(a2 * F0ProOrg[k].at<float>(i,j) + b2);
                    }
                    if ((F0ProOrg[k].at<float>(i,j) >= x2) && (F0ProOrg[k].at<float>(i,j) <= x3)){
                        select2[k].at<int>(i,j) = 1;
                    } else{
                        select2[k].at<int>(i,j) = 0;
                    }
                }
            }
        }
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    if (select2[k].at<int>(i,j) ==1){
                        varMapOut[k].at<float>(i,j) = static_cast<float>(a3 * F0ProOrg[k].at<float>(i,j) + b3);
                    }
                    if (F0ProOrg[k].at<float>(i,j) >= x3){
                        varMapOut[k].at<float>(i,j) = static_cast<float>(y3);
                    }
                }
            }
        }
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    if (isnan(varMapOut[k].at<float>(i,j))){
                        varMapOut[k].at<float>(i,j) = static_cast<float>(y0);
                    }
                }
            }
        }


        return varMapOut;
    }//fit_F0-var


    vector<vector<cv::Mat>> movmean(const vector<vector<cv::Mat>>& dataIn) {
        int H = dataIn[0][0].rows;
        int W = dataIn[0][0].cols;
        int L = dataIn[0].size();
        int T = dataIn.size();
        vector<vector<cv::Mat>> dataOut(T, vector<cv::Mat>(L));
        int halfWin = opts.movAvgWin/2;
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            for (int t = 0; t < T; ++t) {
                int start = max(0, t - halfWin);
                int end = min(T - 1, t + halfWin);
                int count = 0;
                cv::Mat sum = cv::Mat::zeros(dataIn[0][0].rows, dataIn[0][0].cols, CV_32F);
                for (int i = start; i <= end; ++i) {
                    cv::add(sum, dataIn[i][k], sum, cv::noArray(), -1);
                    ++count;
                }
                dataOut[t][k] = sum / count;
            }
        }
        return dataOut;
    }


    vector<cv::Mat> truncated_kept_var(const vector<cv::Mat>& quantiles){
        int H = quantiles[0].rows;
        int W = quantiles[0].cols;
        int L = quantiles.size();
        vector<cv::Mat> pars(L);
        boost::math::normal_distribution<float> normal_dist;
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            pars[k] = cv::Mat(H,W,CV_32F);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    float quantile = quantiles[k].at<float>(i,j);
                    if (quantile == 0.0f){
                        pars[k].at<float>(i,j) = 2.0f;
                    }
                    else{
                        float a = boost::math::quantile(normal_dist, quantile);
                        float phi_a = boost::math::pdf(normal_dist,a);
                        float mu = a * quantile + phi_a;
                        float second_order = a * a * quantile + 1 - quantile + a * phi_a;
                        pars[k].at<float>(i,j) = (second_order - mu * mu) * 2;
                    }
                }
            }
        }//for(k)

        return pars;
    }//truncated_kept_var()

    float truncated_kept_var(float quantile) {
        boost::math::normal_distribution<float> normal_dist;

        if (quantile == 0.0f) {
            return 2.0f;
        }
        else {
            float a = boost::math::quantile(normal_dist, quantile);
            float phi_a = boost::math::pdf(normal_dist, a);
            float mu = a * quantile + phi_a;
            float second_order = a * a * quantile + 1 - quantile + a * phi_a;
            return (second_order - mu * mu) * 2;
        }
    }


    vector<vector<cv::Mat>> baselineLinearEstimate(vector<vector<cv::Mat>>& data){
        int H = data[0][0].rows;
        int W = data[0][0].cols;
        int L = data[0].size();
        int T = data.size();
        vector<vector<cv::Mat>> datMA(T, vector<cv::Mat>(L));
        datMA = movmean(data);
//        for (int t = 0; t < T; t++) {
//            for (int k = 0; k < L; k++) {
//                cv::blur(data[t][k], datMA[t][k], cv::Size(opts.movAvgWin, opts.movAvgWin));
//            }
//        }
//#pragma omp parallel for collapse(4)
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        if (isnan(data[t][k].at<float>(i,j))){
                            datMA[t][k].at<float>(i,j) = NAN;
                        }
                    }
                }
            }
        }
//        cout<<"datMA: "<< endl;
////        for (int k = 0; k < 2; ++k) {
//            for (int i = 0; i < 7; ++i) {
//                for (int j = 0; j < 7; ++j) {
//                    cout<<datMA[0][0].at<float>(i,j)<<" ";
//                }
//                cout<<endl;
//            }
//        }

        int step = static_cast<int>(round(0.5 * opts.cut));
        float maxV = 0;
        double maxVal = 0;
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                cv::minMaxLoc(cv::Mat(data[t][k]), nullptr, &maxVal);
                if (maxVal > maxV){
                    maxV = static_cast<float>(maxVal);
                }
            }
        }
        int nSegment = static_cast<int>(max(1.0, ceil(T/step)-1));
//        vector<vector<cv::Mat>> minPosition(nSegment, vector<cv::Mat>(L, cv::Mat::ones(H, W, CV_32S)));
        vector<vector<cv::Mat>> minPosition(nSegment, vector<cv::Mat>(L));
////#pragma omp parallel for collapse(2)
        for (int kk = 0; kk < nSegment; ++kk) {
            int t0 = kk * step;
            int t1 = min(T, t0 + opts.cut);
            for (int k = 0; k < L; ++k) {
                minPosition[kk][k] = cv::Mat::ones(H, W, CV_32S);
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        float minV = 10;
                        for (int t = t0; t < t1; ++t) {
                            if (datMA[t][k].at<float>(i, j) < minV) {
                                minV = datMA[t][k].at<float>(i, j);
                                minPosition[kk][k].at<int>(i, j) = t;
                            }//if
                        }//for(t)
                    }
                }
            }//for(k)
//            float minV = 10;
//            for (int t = t0; t < t1; ++t) {
//                if (datMA[t][0].at<float>(1,1) < minV) {
//                    minV = datMA[t][0].at<float>(1,1);
//                    minPosition[kk][0].at<int>(1,1) = t;
//                }
//            }
        }//for(kk)

//        for (int k = 0; k < 10; ++k) {
//            cout<<k<<" :"<< minPosition[0][k].at<int>(1,1)<<"  "<< endl;
//        }
//        cout<<"minPosition"<<endl;
//        cout<<minPosition[0][0].at<int>(1,1)<<" : ";
//        cout<<datMA[minPosition[0][0].at<int>(1,1)][0].at<float>(1,1)<<"   ";
//        cout<<endl;
//        cout<<"datMA[t]"<<endl;
//        for (int t = 0; t < 200; ++t) {
//            cout<<t<<" :"<< datMA[t][0].at<float>(1,2)<<"  ";
//        }
//        cout<<"minPosition"<<endl;
//        for (int i = 1; i < 20; ++i) {
//            cout<<minPosition[0][0].at<int>(i,1)<<"  ";
//        }
//        cout<<endl;


        vector<vector<cv::Mat>> F0(T, vector<cv::Mat>(L));
//#pragma omp parallel for collapse(2)
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                F0[t][k] = cv::Mat::ones(H, W, CV_32F);
            }
        }
        for (int k = 0; k < L; ++k) {
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    vector<int> curP;
                    vector<float> value;
                    for (int kk = 0; kk < nSegment; ++kk) {
                        curP.push_back(minPosition[kk][k].at<int>(i,j));
                    }
                    if (nSegment > 1){//try unsorted_set
                        sort(curP.begin(),curP.end());
                        auto it = unique(curP.begin(),curP.end());
                        curP.erase(it, curP.end());
                    }
                    for (int tt = 0; tt < curP.size(); ++tt) {
                        value.push_back(datMA[curP[tt]][k].at<float>(i,j));
//                            if (!isnan(value[tt])){
//                                curP[tt] = curP[tt];
//                                value[tt] = value[tt];
//                            }
                    }

                    int nMin = value.size();
                    cv::Mat curve = cv::Mat::zeros(1,T,CV_32F);
                    if (nMin == 0){
                        curve = maxV;
                    }
                    else{
                        //first part
//#pragma omp parallel for
                        for (int l = 0; l < curP[0]; ++l) {
                            curve.at<float>(0,l) = value[0];
                        }
                        //end part
//#pragma omp parallel for
                        for (int l = curP[nMin-1]; l < T; ++l) {
                            curve.at<float>(0,l) = value[nMin-1];
                        }
                        //middle part
//#pragma omp parallel for
                        for (int l = 0; l < nMin-1; ++l) {
                            int mt1 = curP[l];
                            int mt2 = curP[l+1];
                            for (int m = mt1, index=0; m < mt2; ++m, ++index) {
                                curve.at<float>(0,m) = value[l] + (value[l+1]-value[l])/static_cast<float>((mt2-mt1)*index);
                            }

                        }
                    }//else
//#pragma omp parallel for
                    for (int t = 0; t < T; ++t) {
                        F0[t][k].at<float>(i,j) = curve.at<float>(0,t);
                    }
                }//for(j)
            }//for(i)
        }//for(k)
//#pragma omp parallel for collapse(4)
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        if (isnan(data[t][k].at<float>(i,j))){
                            F0[t][k].at<float>(i,j) = maxV;
                        }
                    }
                }
            }
        }

//        cout<<"F0: "<< endl;
//        cout<<"height of image:"<< F0[0][0].rows << endl;
//        cout<<"width of image:"<< F0[0][0].cols << endl;
//        cout<<"length of image:"<< F0[0].size() << endl;
//        cout<<"time frames of image:"<< F0.size() << endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<F0[0][0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }

        return F0;
    }//baselineLinearEstimate()


    void correctBoundaryStd(){
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

        vector<cv::Mat> correctMap(dist*2+1);
//#pragma omp parallel for
        for (int k = 0; k < dist * 2 + 1; ++k) {
            correctMap[k] = cv::Mat::zeros(dist*2+1, dist*2+1, CV_32F);
        }

//        for (int k = dist; k < 2*dist+1; ++k) {
//            for (int i = dist; i < 2*dist+1; ++i) {
//                for (int j = dist; j < 2*dist+1; ++j) {
//                    vector<cv::Mat> filter1(dist*2+1);
//                    for (int kk = 0; kk < dist * 2 + 1; ++kk) {
//                        filter1[k] = filter0[k];
//                    }//for(kk)
//                    for (int x = ; x < ; ++x) {
//
//                    }
//                }
//            }
//        }

    }


    float obtainBias(){
        MATFile *pmatFile;
        mxArray *pMxbiasMatrix; //215*55 double
        mxArray *pMxcuts; //1*55 double
        mxArray *pMxwindowSizes; //1*215 double
        double *biasMatrix;
        double *cuts;
        double *windowSizes;

//        cout<< "--------loading cfg F0_biasMatrix--------"<<endl;
        const char *filename = "../cfg/F0_biasMatrix.mat";
        pmatFile = matOpen(filename, "r");
        if (pmatFile == nullptr) {
            cout<< "--------error opening cfg--------"<<endl;
            exit(-1);
        }

        pMxbiasMatrix = matGetVariable(pmatFile, "biasMatrix");
        if (pMxbiasMatrix == nullptr) {
            cout<< "--------error reading variable \"biasMatrix\" from file--------"<<endl;
            exit(-1);
        }
        pMxcuts = matGetVariable(pmatFile, "cuts");
        if (pMxcuts == nullptr) {
            cout<< "--------error reading variable \"cuts\" from file--------"<<endl;
            exit(-1);
        }
        pMxwindowSizes = matGetVariable(pmatFile, "windowSizes");
        if (pMxwindowSizes == nullptr) {
            cout<< "--------error reading variable \"windowSizes\" from file--------"<<endl;
            exit(-1);
        }

        biasMatrix = mxGetPr(pMxbiasMatrix);
        if (biasMatrix == nullptr) {
            cout<< "--------error reading data from variable \"biasMatrix\"-------"<<endl;
            exit(-1);
        }
        cuts = mxGetPr(pMxcuts);
        if (cuts == nullptr) {
            cout<< "--------error reading data from variable \"cuts\"-------"<<endl;
            exit(-1);
        }
        windowSizes = mxGetPr(pMxwindowSizes);
        if (windowSizes == nullptr) {
            cout<< "--------error reading data from variable \"windowSizes\"-------"<<endl;
            exit(-1);
        }

        const mwSize *dims_biasMatrix = mxGetDimensions(pMxbiasMatrix);
        const mwSize *dims_cuts = mxGetDimensions(pMxcuts);
        const mwSize *dims_windowSizes = mxGetDimensions(pMxwindowSizes);

        float bias=0;
//        H = dims[0];
//        cout<<"original size: "<< endl;
//        frame[t][k].at<float>(i,j) = static_cast<float>(pdata[j_src*H + i_src + k*H*W + t*H*W*L]);
        if (opts.movAvgWin > opts.cut){
            bias = 0;
            return bias;
        }

        //linear interpolate value if cannot find it
        int idx0=-1, idy0=-1, idx1=-1, idy1=-1;
//#pragma omp parallel for
        for (int i = dims_windowSizes[1] - 1; i >= 0; --i) {
            if (windowSizes[i] <= opts.movAvgWin){
                idx0 = i;
                break;
            }
        }
//#pragma omp parallel for
        for (int i = dims_cuts[1] - 1; i >= 0; --i) {
            if (cuts[i] <= opts.cut){
                idy0 = i;
                break;
            }
        }
        if ((idx0==dims_windowSizes[1]-1) || (windowSizes[idx0]==opts.movAvgWin)){
            idx1 = idx0;
        } else{
            idx1 = idx0 + 1;
        }
        if ((idy0==dims_cuts[1]-1) || (cuts[idy0]==opts.cut)){
            idy1 = idy0;
        } else{
            idy1 = idy0 + 1;
        }

        float bias0=0, bias1=0;
//#pragma omp parallel for
        for (int i = idx0; i <= idx1 ; ++i) {
            if (!isnan(biasMatrix[i + idy0 * dims_biasMatrix[0]])){
                bias0 += static_cast<float>(biasMatrix[i + idy0 * dims_biasMatrix[0]]);
            }
            if (!isnan(biasMatrix[i + idy1 * dims_biasMatrix[0]])){
                bias1 += static_cast<float>(biasMatrix[i + idy1 * dims_biasMatrix[0]]);
            }
        }
        bias0 /= static_cast<float>(idx1 - idx0 + 1);
        bias1 /= static_cast<float>(idx1 - idx0 + 1);
        int cut0 = static_cast<int>(cuts[idy0]);
        int cut1 = static_cast<int>(cuts[idy0+1]);

        if (isnan(bias0)){
            bias = bias1;
        } else{
            bias = bias0 + (bias1 - bias0) / (cut1 - cut0) * (opts.cut - cut0);
        }

        return bias;
    }//obtainBias()


    void noiseEstimationFunction(const vector<vector<cv::Mat>>& dataOrg, const vector<vector<cv::Mat>>& dataSmo,
                                 const vector<cv::Mat>& F0Pro, bool*** evtSpatialMask, vector<cv::Mat>& stdMapOrg, vector<cv::Mat>& stdMapSmo,
                                 vector<cv::Mat>& tempVarOrg,   vector<cv::Mat>& correctPars){
        /*
         * return "stdMapOrg,stdMapSmo,tempVarOrg,correctPars"
         * variance map
         * calculate the variance of raw data
         */
        int H = dataOrg[0][0].rows;
        int W = dataOrg[0][0].cols;
        int L = dataOrg[0].size();
        int T = dataOrg.size();
        bool correctNoise = true;
        vector<cv::Mat> tempMap(L);
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            tempMap[k] = cv::Mat(H,W,CV_32F);
            tempVarOrg[k] = cv::Mat(H,W,CV_32F);
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    double sum = 0;
                    for (int t = 0; t < T-1; ++t) {
                        sum += pow((dataOrg[t][k].at<float>(i,j) - dataOrg[t+1][k].at<float>(i,j)), 2);
                    }
                    tempMap[k].at<float>(i,j) = static_cast<float>(sum/T);
                }
            }
            tempVarOrg[k] = tempMap[k] / 2;
        }//for(k)
//        cout<<"tempMap: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<tempMap[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
        vector<cv::Mat> varMapOrg(L);
        if (correctNoise){
            vector<cv::Mat> countInValid(L);
            vector<cv::Mat> totalSamples(L);
            vector<cv::Mat> ratio(L);
//#pragma omp parallel for collapse(3)
            for (int k = 0; k < L; ++k) {
                countInValid[k] = cv::Mat(H,W,CV_32F);
                totalSamples[k] = cv::Mat(H,W,CV_32F);
                ratio[k] = cv::Mat(H,W,CV_32F);
                for (int i = 0; i < H; ++i) {
                    for (int j = 0; j < W; ++j) {
                        int count_invalid = 0, count_samples=0;
                        for (int t = 0; t < T; ++t) {
                            if (dataOrg[t][k].at<float>(i,j) == 0){
                                ++count_invalid;
                            }
                            if (!isnan(dataOrg[t][k].at<float>(i,j))){
                                ++count_samples;
                            }
                            countInValid[k].at<int>(i,j) = count_invalid;
                            totalSamples[k].at<int>(i,j) = count_samples;
                        }
                    }
                }
                ratio[k] = countInValid[k] / totalSamples[k];
            }//for(k)
            correctPars = truncated_kept_var(ratio);
//#pragma omp parallel for
            for (int k = 0; k < L; ++k) {
                varMapOrg[k] = cv::Mat(H,W,CV_32F);
                varMapOrg[k] = tempMap[k] / correctPars[k];
            }
        }//if(correctNoise)
        else{
//#pragma omp parallel for
            for (int k = 0; k < L; ++k) {
                varMapOrg[k] = cv::Mat(H,W,CV_32F);
                varMapOrg[k] = tempVarOrg[k];
            }
        }//else
//#pragma omp parallel for collapse(3)
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    if (!evtSpatialMask[i][j][k]){
                        varMapOrg[k].at<float>(i,j) = NAN;
                    }
                }
            }
        }
//        cout<<"F0Pro: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<F0Pro[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
//        cout<<"varMapOrg: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<varMapOrg[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
        vector<cv::Mat> varMapOut = fit_F0_var(F0Pro, varMapOrg, H, W, L, T);
//        cout<<"varMapOut: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<varMapOut[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
        vector<cv::Mat> varMapSmo(L);
//#pragma omp parallel for
        for (int k = 0; k < L; ++k) {
            cv::sqrt(varMapOut[k],stdMapOrg[k]);
            varMapSmo[k] = cv::Mat(H,W,CV_32F);
        }
        if (opts.smoXY == 0){
//#pragma omp parallel for
            for (int k = 0; k < L; ++k) {
                stdMapSmo[k] = stdMapOrg[k];
                varMapSmo[k] = varMapOrg[k];
            }
        }
//        cout<<"stdMapOrg: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<stdMapOrg[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
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
        //estimated variance from smoothed data
//#pragma omp parallel for collapse(3)
        for (int k = 0; k < L; ++k) {
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    double sum = 0;
                    for (int t = 0; t < T-1; ++t) {
                        sum += pow(dataSmo[t][k].at<float>(i,j) - dataSmo[t+1][k].at<float>(i,j), 2);
                    }
                    varMapSmo[k].at<float>(i,j) = static_cast<float>(sum/T);
                }
            }
        }//for(k)
//        cout<<"varMapSmo: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<varMapSmo[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
        //correct the variance according to truncated model
        if (correctNoise){
//#pragma omp parallel for
            for (int k = 0; k < L; ++k) {
                cv::Mat temp1 = cv::Mat(H,W,CV_32F);
                cv::Mat temp2 = cv::Mat(H,W,CV_32F);
                cv::filter2D(tempMap[k]/correctPars[k], temp1, -1, filter[dist]);
                cv::filter2D(varMapOrg[k], temp2, -1, filter[dist]);
                cv::divide(temp1, temp2, temp1);
                cv::multiply(varMapSmo[k], temp1, varMapSmo[k]);
//                varMapSmo[k] = varMapSmo[k] * temp1 / temp2;
            }
        }//if
//#pragma omp parallel for collapse(3)
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    if (!evtSpatialMask[i][j][k]){
                        varMapSmo[k].at<float>(i,j) = NAN;
                    }
                }
            }
        }
        //correct the variance in the boundary(caused by smoothing operation)
//        correctMap2 = correctBoundaryStd();
//        vector<cv::Mat> correctMap2(L);
//        for (int k = 0; k < L; ++k) {
//            correctMap2[k] = cv::Mat::zeros(H,W,CV_32S);
//            for (int i = dist; i < H-dist; ++i) {
//                for (int j = dist; j < W-dist; ++j) {
//                    correctMap2[k].at<int>(i,j) = 1;
//                }
//            }
//        }
//        for (int k = 0; k < L; ++k) {
//            for (int i = 0; i < H; ++i) {
//                for (int j = 0; j < W; ++j) {
//                    if (correctMap2[k].at<int>(i,j) == 0){
//                        varMapSmo[k].at<float>(i,j) = NAN;
//                    }
//                }
//            }
//        }
//        writeDataToMatFile(F0Pro,"C:/Users/Kevin Qiao/Desktop/AQuA_data/F0pro.mat");
//        writeDataToMatFile(varMapSmo,"C:/Users/Kevin Qiao/Desktop/AQuA_data/varMapSmo.mat");
//        cout<<"finish";
//        cout<<"varMapSmo: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<varMapSmo[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
        varMapOut = fit_F0_var(F0Pro, varMapSmo, H, W, L, T, dist);
//#pragma omp parallel for
        for (int k = 0; k < L; ++k) {
            cv::sqrt(varMapOut[k],stdMapSmo[k]);
        }

    }//noiseEstimationFunction()


    vector<vector<cv::Mat>> baselineRemoveAndNoiseEstimation(vector<vector<cv::Mat>>& dataOrg, bool*** evtSpatialMask){
        /*
         * smooth the data
         */
        int H = dataOrg[0][0].rows;
        int W = dataOrg[0][0].cols;
        int L = dataOrg[0].size();
        int T = dataOrg.size();
        cout<< "--------start baselineRemoveAndNoiseEstimation--------"<<endl;
        int ksize = 2 * ceil(2 * opts.smoXY) + 1;
        vector<vector<cv::Mat>> dataSmo(T, vector<cv::Mat>(L));
        if (opts.smoXY > 0 ){
//#pragma omp parallel for collapse(2)
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    cv::GaussianBlur(dataOrg[t][k],dataSmo[t][k],cv::Size(ksize, ksize), opts.smoXY, opts.smoXY);
                }
//                // 3rd dimension smoothing
//                for (int i = 0; i < H; ++i) {
//                    for (int j = 0; j < W; ++j) {
//                        cv::Mat slice(L,1,CV_32F);
//                        for (int k = 0; k < L; ++k) {
//                            slice.at<float>(k) = dataSmo[t][k].at<float>(i,j);
//                        }
//                        cv::GaussianBlur(slice,slice,cv::Size(1,ksize),0,opts.smoXY);
//                        for (int k = 0; k < L; ++k) {
//                            dataSmo[t][k].at<float>(i,j) = slice.at<float>(k);
//                        }
//                    }
//                }
            }//for(t)
        }//if(smoXY>0)

//        cout<<"dataSmo: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<dataSmo[0][0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }


        /*
         * linear estimation of F0
         */
        opts.cut = min(opts.cut, T);
        //remove baseline
        vector<vector<cv::Mat>> F0(T, vector<cv::Mat>(L));
        vector<cv::Mat> F0Pro(L);
//#pragma omp parallel for
        for (int k = 0; k < L; ++k) {
            F0Pro[k] = cv::Mat(H,W,CV_32F);
        }
        F0 = baselineLinearEstimate(dataSmo);
//#pragma omp parallel for collapse(3)
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    float sum = 0;
                    for (int t = 0; t < T; ++t) {
                        sum += F0[t][k].at<float>(i,j);
                    }
                    F0Pro[k].at<float>(i,j) = sum/T;
                }
            }
        }

        //noise estimation
        vector<cv::Mat> stdMapOrg(L);
        vector<cv::Mat> stdMapGau(L);
        vector<cv::Mat> tempVarOrg(L);
        vector<cv::Mat> correctPars(L);
        noiseEstimationFunction(dataOrg, dataSmo, F0Pro, evtSpatialMask, stdMapOrg, stdMapGau, tempVarOrg, correctPars);
//        cout<<"stdMapOrg: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<stdMapOrg[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
//        cout<<"stdMapGau: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<stdMapGau[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
//        cout<<"tempVarOrg: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<tempVarOrg[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
//        cout<<"correctPars: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<correctPars[0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
//        cout<<"F0: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<F0[0][0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
//        cout<<"dataSmo: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<dataSmo[0][0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
        float bias = obtainBias();
//        cout<<"bias: "<< bias<<endl;
        vector<vector<cv::Mat>> dF(T, vector<cv::Mat>(L));
        // correct bias during noise estimation. Bias does not impact noise
        // normalization - zScoreMap
//#pragma omp parallel for collapse(2)
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                dF[t][k] = cv::Mat(H, W, CV_32F);
                cv::subtract(F0[t][k], stdMapGau[k] * bias, F0[t][k]);
                cv::subtract(dataSmo[t][k], F0[t][k], dF[t][k]);
                cv::divide(dF[t][k], stdMapGau[k], dF[t][k]);
            }
        }
//        cout<<"dF: "<< endl;
//        for (int i = 0; i < 7; ++i) {
//            for (int j = 0; j < 7; ++j) {
//                cout<<dF[0][0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }

        //if (ch==1)
        opts.stdMapGau1 = stdMapGau;
        opts.stdMapOrg1 = stdMapOrg;
        opts.tempVarOrg1 = tempVarOrg;
        opts.correctPars1 = correctPars;
//        for (int k = 0; k < L; ++k) {
//            opts.stdMapGau1.push_back(stdMapGau[k].clone());
//            opts.stdMapOrg1.push_back(stdMapOrg[k].clone());
//            opts.tempVarOrg1.push_back(tempVarOrg[k].clone());
//            opts.correctPars1.push_back(correctPars[k].clone());
//        }
//        cout<< opts.stdMapOrg1.size()<<endl;
//        cout<< opts.stdMapOrg1[0].size()<<endl;

        //for visualize
        double maxVal = 0;
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                cv::minMaxLoc(cv::Mat(dF[t][k]), nullptr, &maxVal);
                if (maxVal > opts.maxdF1){
                    opts.maxdF1 = static_cast<float>(maxVal);
                }
            }
        }

        return dF;
    }//baselineRemoveAndNoiseEstimation()


}// namespace