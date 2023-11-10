//
// Created by Kevin Qiao on 7/18/2023.
//

#include "phaseRun.h"

namespace AQuA{

    cv::Mat myResize(const cv::Mat& src, float scaleRatio_x, float scaleRatio_y){
        float h = ceil(src.rows/scaleRatio_y); //height of dst
        float w = ceil(src.cols/scaleRatio_x); // width of dst
        int unit_h = ceil(src.rows / h); //number of units in row
        int unit_w = ceil(src.cols / w); // number of units in col
        cv::Mat dst = cv::Mat(h, w, CV_32F);
        for (int i_dst = 0; i_dst < h; ++i_dst) {
            for (int j_dst = 0; j_dst < w; ++j_dst) {
                float sum = 0;
                float count = 0;
                for (int i_src = i_dst*unit_h; (i_src<(i_dst+1)*unit_h) && (i_src<src.rows); ++i_src) {
                    for (int j_src = j_dst*unit_w; (j_src<(j_dst+1)*unit_w) && (j_src<src.cols); ++j_src) {
                        sum += src.at<float>(i_src,j_src);
                        ++count;
//                        cout<<"i_src: "<<i_src<<"  j_src: "<<j_src<<" ";
                    }//for(j_src)
//                    cout<<endl;
                }//for(i_src)
                dst.at<float>(i_dst, j_dst) = sum/count;
            }
        }
        return dst;
    }//myResize()


    void ordStatSmallSampleWith0s(const vector<float>& fg, const vector<float>& bg, const vector<float>& nanVec, double& mu, double& sigma){
        static vector<vector<double>> mus;
        static vector<cv::Mat> covMatrixs;

        if(fg.empty() && bg.empty()){
            mu = NAN;
            sigma = NAN;
            return;
        }

        int M = fg.size();
        int N = bg.size();
        int nanLen = nanVec.size();
        int n = M + N + nanLen;

        float** all = new float* [n];
        for (int i = 0; i < n; ++i) {
            all[i] = new float [2];    //all[i][0] = all; all[i][1] = label;
        }

        for (int i = 0; i < N; ++i) {
            all[i][0] = bg[i];
            all[i][1] = -1;
        }
        for (int i = N, j=0; i < M+N; ++i,++j) {
            all[i][0] = fg[j];
            all[i][1] = 1;
        }
        for (int i = M+N, j=0; i < n; ++i,++j) {
            all[i][0] = nanVec[j];
            all[i][1] = 0;
        }

        auto compareFunc = [](const float* a, const float* b) {
            return a[0] < b[0];
        };
        sort(all,all+n, compareFunc);

        if (mus.empty() || covMatrixs.empty()){
            mus = loadCell_double("../cfg/Order_mus_sigmas.mat","mus");
            covMatrixs = loadCell_matrix("../cfg/Order_mus_sigmas.mat","covMatrixs");
        }

        vector<double> muVec = mus[n-1];
        vector<int> ind1;
        vector<int> indm1;
        vector<int> ind0;
        for (int i = 0; i < muVec.size(); ++i) {
            if (all[i][1] == 1){
                ind1.emplace_back(i);
            }
            if (all[i][1] == -1){
                indm1.emplace_back(i);
            }
            if (all[i][1] == 0){
                ind0.emplace_back(i);
            }
        }
        double tempSum1 = 0;
        double tempSum2 = 0;
        int cnt1 = 0;
        int cnt2 = 0;
        for (int inx: ind1) {
            tempSum1+=muVec[inx];
            ++cnt1;
        }
        for (int inx: indm1) {
            tempSum2+=muVec[inx];
            ++cnt2;
        }
        mu = tempSum1/cnt1 - tempSum2/cnt2;
        cv::Mat covMatrix = covMatrixs[n-1].clone();
        for (int inx: ind0) {
            covMatrix.row(inx) = 0;
            covMatrix.col(inx) = 0;
        }
        for (int inx: ind1) {
            cv::divide(covMatrix.row(inx),static_cast<float>(M),covMatrix.row(inx));
            cv::divide(covMatrix.col(inx),static_cast<float>(M),covMatrix.col(inx));
        }
        for (int inx: indm1) {
            cv::divide(covMatrix.row(inx),static_cast<float>(-N),covMatrix.row(inx));
            cv::divide(covMatrix.col(inx),static_cast<float>(-N),covMatrix.col(inx));
        }

        sigma = sqrt(static_cast<float>(cv::sum(covMatrix)[0]));
        return;
    }


    Score_struct getSeedScore_DS4(const vector<int>& pix, const vector<vector<cv::Mat>>& datVec, int H, int W, int L, int T, float t_scl){
        Score_struct result;
        //down sample
        vector<int> ih;
        vector<int> iw;
        vector<int> il;
        vector<int> it;
        for (int ii_ite = 0; ii_ite < pix.size(); ++ii_ite) {
            ih.emplace_back(ind2sub(pix[ii_ite],H,W,L).i);
            iw.emplace_back(ind2sub(pix[ii_ite],H,W,L).j);
            il.emplace_back(ind2sub(pix[ii_ite],H,W,L).k);
            it.emplace_back(ind2sub(pix[ii_ite],H,W,L).t);
        }

        for (int i = 0; i < pix.size(); ++i) {
            it[i] = floor(it[i]/t_scl);
        }
        int T0 = floor(T/t_scl);

        cv::Mat select = cv::Mat::ones(pix.size(),1,CV_8U);
        for (int i = 0; i < pix.size(); ++i) {
            if (it[i] >= T0){
                select.at<uchar>(i,0) = 0;
            }
        }

        vector<int> pix0;
        for (int i = 0; i < ih.size(); ++i) {
            if (select.at<uchar>(i,0) == 1) {
                pix0.emplace_back(sub2ind(ih[i],iw[i],il[i],it[i],H,W,L));
            }
        }

        unordered_map<int,int> a_counts;
        for (const int& num : pix0) {
            ++a_counts[num];
        }
        std::vector<int> fgPix;
        for (const auto& pair : a_counts) {
            if (pair.second > t_scl / 2) {
                fgPix.emplace_back(pair.first);
            }
        }
        sort(fgPix.begin(), fgPix.end());

        if (fgPix.empty()){
            result.z_score1 = 0;
            result.z_score2 = 0;
            result.t_score1 = 0;
            result.t_score2 = 0;
            return result;
        }

        //find neighbor
        ih.clear();
        iw.clear();
        il.clear();
        it.clear();
        for (int i = 0; i < fgPix.size(); ++i) {
            ih.emplace_back(ind2sub(fgPix[i],H,W,L).i);
            iw.emplace_back(ind2sub(fgPix[i],H,W,L).j);
            il.emplace_back(ind2sub(fgPix[i],H,W,L).k);
            it.emplace_back(ind2sub(fgPix[i],H,W,L).t);
        }

        vector<int> ihwOrg;
        unordered_set<int> ihw_temp;

        for (int i = 0; i < ih.size(); ++i) {
            ihwOrg.emplace_back(sub2ind(ih[i],iw[i],il[i],H,W));
            ihw_temp.insert(ihwOrg[i]);
        }
        vector<int> ihw(ihw_temp.begin(), ihw_temp.end());
        sort(ihw.begin(), ihw.end());

        vector<int> sz;
        vector<vector<float>> fgAll(ihw.size());
        vector<vector<float>> bgL(ihw.size());
        vector<vector<float>> bgR(ihw.size());
        vector<vector<float>> nanVec(ihw.size());
        int cnt = 0;
        int cnt2 = 0;
        vector<float> noise;
        vector<int> degreeOfFreedoms;

        for (int ii_ihw = 0; ii_ihw < ihw.size(); ++ii_ihw) { // --i
            vector<int> curIt;
            vector<int>::iterator iter = ihwOrg.begin();
            while ((iter = find(iter, ihwOrg.end(), ihw[ii_ihw])) != ihwOrg.end()) {
                curIt.push_back(it[distance(ihwOrg.begin(), iter)]);
                ++iter;
            }


            //get normalized down sampled curve
            cv::Mat curve0 = cv::Mat(T,1,CV_32F);
            int h_datVec = ind2sub(ihw[ii_ihw],datVec[0][0].rows,datVec[0][0].cols).i;
            int w_datVec = ind2sub(ihw[ii_ihw],datVec[0][0].rows,datVec[0][0].cols).j;
            int k_datVec = ind2sub(ihw[ii_ihw],datVec[0][0].rows,datVec[0][0].cols).k;
            #pragma omp parallel for
            for (int t = 0; t < T; ++t) {
                curve0.at<float>(t,0) = datVec[t][k_datVec].at<float>(h_datVec,w_datVec);
            }
            cv::Mat curve = myResize(curve0,t_scl,t_scl);
            cv::multiply(curve,sqrt(t_scl),curve);
            int dur = curIt.size();
            sz.emplace_back(dur);
            int t0 = *min_element(curIt.begin(),curIt.end());
            int t1 = *max_element(curIt.begin(),curIt.end());
            vector<int> t_Left;
            vector<int> t_Right;
            for (int tt = max(0,t0-dur); tt <= t0-1; ++tt) {
                t_Left.emplace_back(tt);
            }
            for (int tt = t1+1; tt <= min(T0, t1 + dur); ++tt) {
                t_Right.emplace_back(tt);
            }
            for (int ite = 0;ite<curIt.size();++ite) {
                fgAll[ii_ihw].emplace_back(curve.at<float>(curIt[ite],0));
            }
            for (int ite = 0;ite<t_Left.size();++ite) {
                bgL[ii_ihw].emplace_back(curve.at<float>(t_Left[ite],0));
            }
            for (int ite = 0;ite<t_Right.size();++ite) {
                bgR[ii_ihw].emplace_back(curve.at<float>(t_Right[ite],0));
            }
            vector<int> diff;
            vector<int> t0_t1;
            for (int ii = t0; ii <= t1; ++ii) {
                t0_t1.emplace_back(ii);
            }
            set_difference(t0_t1.begin(),t0_t1.end(),curIt.begin(),curIt.end(),
                                inserter(diff,diff.begin()));
            for (int ite = 0;ite<diff.size();++ite) {
                nanVec[ii_ihw].emplace_back(curve.at<float>(diff[ite],0));
            }

            //jump one point
            vector<float> difL;
            vector<float> difR;
            if (!t_Left.empty()){
                for (auto& ite: t_Left) {
                    ite -= 1;
                }
                t0 = max(static_cast<float>(1),*min_element(t_Left.begin(),t_Left.end())*t_scl);
                t1 = *max_element(t_Left.begin(),t_Left.end())*t_scl;
                for (int tt = t0;tt <= t1;++tt) {
                    difL.emplace_back(pow((curve0.at<float>(tt,0)-curve0.at<float>(tt-1,0)),2));
                }
                for (int tt = t0-1; tt <= t1; ++tt) {
                    if(curve0.at<float>(tt,0) == 0){
                        ++cnt;
                    }
                }
                cnt2 = cnt2 + t1 -t0 + 2;
            }//if(!t_left.empty())

            if (!t_Right.empty()){
                for (auto& ite: t_Right) {
                    ite += 1;
                }
                t0 = *min_element(t_Right.begin(),t_Right.end())*t_scl;
                t1 = min(static_cast<float>(T-2), *max_element(t_Right.begin(),t_Right.end())*t_scl);
                for (int tt = t0;tt <= t1;++tt) {
                    difR.emplace_back(pow((curve0.at<float>(tt,0)-curve0.at<float>(tt+1,0)),2));
                }
                for (int tt = t0; tt <= t1+1; ++tt) {
                    if(curve0.at<float>(tt,0) == 0){
                        ++cnt;
                    }
                }
                cnt2 = cnt2 + t1 -t0 + 2;
            }//if(!t_Right.empty())

            noise.insert(noise.end(),difL.begin(),difL.end());
            noise.insert(noise.end(),difR.begin(),difR.end());
            degreeOfFreedoms.emplace_back(difL.size());
            degreeOfFreedoms.emplace_back(difR.size());

        }//for(ii_ihw)  --i



        vector<float> fg0;
        vector<float> bkL;
        vector<float> bkR;
        for (int i = 0; i < fgAll.size(); ++i) {
            for (int j = 0; j < fgAll[i].size(); ++j) {
                if (!isnan(fgAll[i][j])){
                    fg0.emplace_back(fgAll[i][j]);
                }
                if (!isnan(bgL[i][j])){
                    bkL.emplace_back(bgL[i][j]);
                }
                if (!isnan(bgR[i][j])){
                    bkR.emplace_back(bgR[i][j]);
                }
            }
        }
        int n1 = bkL.size();
        int n2 = bkR.size();

        if ((n1+n2<=2) || (n1<=1) || (n2<=1) || (noise.size()<=2)){
            result.z_score1 = 0;
            result.z_score2 = 0;
            result.t_score1 = 0;
            result.t_score2 = 0;
            return result;
        }

        float correctPar = truncated_kept_var(cnt/cnt2);
        float noise_sum = 0;
        for (int i = 0; i < noise.size(); ++i) {
            noise_sum += noise[i];
        }
        float sigma0 = sqrt(noise_sum/noise.size()/correctPar)/ sqrt(t_scl);

        vector<double> mus_Left(ihw.size());
        vector<double> L_left(ihw.size());
        vector<double> mus_Right(ihw.size());
        vector<double> L_right(ihw.size());
        for (int i = 0; i < ihw.size(); ++i) {
            mus_Left[i] = 0;
            mus_Right[i] = 0;
            L_left[i] = -INFINITY;
            L_right[i] = -INFINITY;
        }

        for (int i = 0; i < ihw.size(); ++i) {  // for i = 1:numel(ihw)
            vector<float> fg;
            vector<float> bg1;
            vector<float> bg2;
            vector<float> nanV;
            double mu;
            double sigma;
            if (!fgAll[i].empty()){
                for (int ii = 0; ii < fgAll[i].size(); ++ii) {
                    fg.emplace_back(fgAll[i][ii]/sigma0);
                    bg1.emplace_back(bgL[i][ii]/sigma0);
                    bg2.emplace_back(bgR[i][ii]/sigma0);
                }
            }
            if (!nanVec[i].empty()){
                for (int ii = 0; ii < nanVec[i].size(); ++ii) {
                    nanV.emplace_back(nanVec[i][ii]/sigma0);
                }
            }
//            for (int ii = 0; ii < fgAll[i].size(); ++ii) { // access each element in fgAll[i]  ---ii
//                if (!fgAll[i].empty()){
//                    fg.emplace_back(fgAll[i][ii]/sigma0);
//                }
//                if (!bgL[i].empty()){
//                    bg1.emplace_back(bgL[i][ii]/sigma0);
//                }
//                if (!bgR[i].empty()){
//                    bg2.emplace_back(bgR[i][ii]/sigma0);
//                }
//                if (!nanVec[i].empty()){
//                    nanV.emplace_back(nanVec[i][ii]/sigma0);
//                }
//            }
//            if (!nanVec[i].empty()){
//                for (int ii = 0; ii < nanVec.size(); ++ii) {
//                    nanV.emplace_back(nanVec[i][ii]/sigma0);
//                }
//            }
            if (!bg1.empty()){
                float LL;
                vector<float> bg_nanV (bg2.begin(), bg2.end());
                if (fg.size()==1){
                    LL = fg[0] - bg1[0];
                } else{
                    float fg_sum = 0;
                    float bg1_sum = 0;
                    for (int j = 0; j < fg.size(); ++j) {
                        fg_sum += fg[j];
                        bg1_sum += bg1[j];
                    }
                    LL = fg_sum/fg.size() - bg1_sum/fg.size();
                } //else
                bg_nanV.insert(bg_nanV.end(),nanV.begin(), nanV.end());
//                cout<< i <<endl;
                ordStatSmallSampleWith0s(fg,bg1,bg_nanV,mu,sigma);
                mus_Left[i] = mu/sigma;
                L_left[i] = LL/sigma;
            }// if (!bg1.empty())
            if (!bg2.empty()){
                float LL;
                vector<float> bg_nanV (bg1.begin(), bg1.end());
                if (fg.size()==1){
                    LL = fg[0] - bg2[0];
                } else{
                    float fg_sum = 0;
                    float bg2_sum = 0;
                    for (int j = 0; j < fg.size(); ++j) {
                        fg_sum += fg[j];
                        bg2_sum += bg2[j];
                    }
                    LL = fg_sum/fg.size() - bg2_sum/fg.size();
                } //else
                bg_nanV.insert(bg_nanV.end(),nanV.begin(), nanV.end());
                ordStatSmallSampleWith0s(fg,bg2,bg_nanV,mu,sigma);
                mus_Right[i] = mu/sigma;
                L_right[i] = LL/sigma;
            }// if (!bg2.empty())

        }// for i = 1:numel(ihw)
        float degreeOfFreedom_sum = 0;
        for (auto i: degreeOfFreedoms) {
            degreeOfFreedom_sum += i;
        }
        float degreeOfFreedom_val = 2 * pow(degreeOfFreedom_sum,2) / (3 * degreeOfFreedom_sum - degreeOfFreedoms.size());
//        save_vector(L_left,"C:\\Users\\Kevin Qiao\\Desktop\\AQuA_data\\test\\L_left.bin");
//        save_vector(L_right,"C:\\Users\\Kevin Qiao\\Desktop\\AQuA_data\\test\\L_right.bin");
//        save_vector(mus_Left,"C:\\Users\\Kevin Qiao\\Desktop\\AQuA_data\\test\\mus_Left.bin");
//        save_vector(mus_Right,"C:\\Users\\Kevin Qiao\\Desktop\\AQuA_data\\test\\mus_Right.bin");
        vector<double>z_Left(L_left.size());
        vector<double>z_Right(L_right.size());
        //  #### non-central cdf is different from the one in matlab
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

        float sz_sum = 0;
        for (int i = 0; i < sz.size(); ++i) {
            result.z_score1 += sqrt(sz[i])*z_Left[i];
            result.z_score2 += sqrt(sz[i])*z_Right[i];
            sz_sum += sz[i];
        }
        result.z_score1 /= sqrt(sz_sum);
        result.z_score2 /= sqrt(sz_sum);

        result.t_score1 = 0;
        for (int i = 0; i < fg0.size(); ++i) {
            fg0[i] /= sigma0;
            bkL[i] /= sigma0;
            bkR[i] /= sigma0;
        }
        if (!bkL.empty()){
            float fg0_sum = 0;
            float bkL_sum = 0;
            for (int i = 0; i < fg0.size(); ++i) {
                fg0_sum += fg0[i];
                bkL_sum += bkL[i];
            }
            result.t_score1 = (fg0_sum/fg0.size() - bkL_sum/fg0.size()) / sqrt(1.0/fg0.size()+1.0/bkL.size());
        }

        result.t_score2 = 0;
        if (!bkR.empty()){
            float fg0_sum = 0;
            float bkR_sum = 0;
            for (int i = 0; i < fg0.size(); ++i) {
                fg0_sum += fg0[i];
                bkR_sum += bkR[i];
            }
            result.t_score2 = (fg0_sum/fg0.size() - bkR_sum/fg0.size()) / sqrt(1.0/fg0.size()+1.0/bkR.size());
        }

        if (noise.size()<100){
            double dof = static_cast<double>(noise.size()) - 1.0;
            boost::math::students_t_distribution<double> student_t(dof);
            double upperCDF1 = 1.0 - boost::math::cdf(student_t, result.t_score1);
            double upperCDF2 = 1.0 - boost::math::cdf(student_t, result.t_score2);
            if (upperCDF1 > 0.999999999999) {
                upperCDF1 = 0.999999999999;
            } else if (upperCDF1 < 0.000000000001) {
                upperCDF1 = 0.000000000001;
            }
            if (upperCDF2 > 0.999999999999) {
                upperCDF2 = 0.999999999999;
            } else if (upperCDF2 < 0.000000000001) {
                upperCDF2 = 0.000000000001;
            }
            boost::math::normal_distribution<double> norm_noise;
            result.t_score1 = -boost::math::quantile(norm_noise, upperCDF1);
            result.t_score2 = -boost::math::quantile(norm_noise, upperCDF2);
        }

//        cout<<"result: "<<result.z_score1<<" "<<result.z_score2;
        return result;
    }//getSeedScore_DS4



    vector<vector<vector<cv::Mat>>> normalizeAndResize(const vector<vector<cv::Mat>>& dataOrg){
        int H = dataOrg[0][0].rows;
        int W = dataOrg[0][0].cols;
        int L = dataOrg[0].size();
        int T = dataOrg.size();
        vector<float> scaleRatios = opts.scaleRatios;
        vector<vector<cv::Mat>> datDS(T,vector<cv::Mat>(L));
        vector<vector<vector<cv::Mat>>> datResize(scaleRatios.size(), vector<vector<cv::Mat>>(T, vector<cv::Mat>(L)));

        //reScl dFOrg
        for (int ii = 0; ii < scaleRatios.size(); ++ii) {
            float scaleRatio = scaleRatios[ii];
//            #pragma omp parallel for collapse(2)
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
//                    cv::resize(dataOrg[t][k],datDS[t][k],cv::Size(),1/scaleRatio,1/scaleRatio, cv::INTER_AREA);
                    datDS[t][k] = myResize(dataOrg[t][k],scaleRatio,scaleRatio);
                }//for(k)
            }//for(t)

            //consider the possible noise correlation, need to re-estimate noise
            vector<cv::Mat> curVarMap(L);
//            #pragma omp parallel for
            for (int k = 0; k < L; ++k) {
                curVarMap[k] = cv::Mat(cv::Mat(datDS[0][0].rows,datDS[0][0].cols,CV_32F));
                cv::Mat sum = cv::Mat::zeros(datDS[0][0].rows,datDS[0][0].cols,CV_32F);
                cv::Mat temp = cv::Mat::zeros(datDS[0][0].rows,datDS[0][0].cols,CV_32F);
                for (int t = 1; t < T; ++t) {
                    cv::subtract(datDS[t][k],datDS[t-1][k],temp);
                    cv::pow(temp,2,temp);
                    cv::add(temp,sum,sum);
                }
                cv::divide(sum,2*(T-1),curVarMap[k]);
//                for (int i = 0; i < datDS[0][0].rows; ++i) {
//                    for (int j = 0; j < datDS[0][0].cols; ++j) {
//                        double sum = 0;
//                        for (int t = 1; t < T; ++t) {
//                            sum += pow((datDS[t][k].at<float>(i,j) - datDS[t-1][k].at<float>(i,j)), 2);
//                        }
//                        curVarMap[k].at<float>(i,j) = static_cast<float>(sum/T/2);
//                    }
//                }
            }//for(k)

            //correct noise due to truncation
            vector<cv::Mat> var1(L);
            vector<cv::Mat> var2(L);
            vector<cv::Mat> curStdMap(L);
//            #pragma omp parallel for
            for (int k = 0; k < L; ++k) {
                var1[k] = myResize(opts.tempVarOrg1[k],scaleRatio,scaleRatio);

                cv::Mat temp1;
                cv::multiply(opts.tempVarOrg1[k], 2, temp1);
                cv::divide(temp1, opts.correctPars1[k], temp1);
//                cv::resize(opts.tempVarOrg1[k], var1[k], cv::Size(), 1 / scaleRatio, 1 / scaleRatio, cv::INTER_AREA);
//                cv::resize(temp1, var2[k], cv::Size(), 1 / scaleRatio, 1 / scaleRatio, cv::INTER_AREA);
                var2[k] = myResize(temp1, scaleRatio,scaleRatio);

                cv::Mat temp2;
                cv::multiply(curVarMap[k], var2[k], temp2);
                cv::divide(temp2, var1[k], temp2);
                cv::sqrt(temp2, curStdMap[k]);
                for (int t = 0; t < T; ++t) {
                    cv::divide(datDS[t][k], curStdMap[k], datDS[t][k]); // since it is only used in seed detection,
                    // and seed detection only check single pixel's score. Not need to fitting again
                }
            }//for(k)

//            #pragma omp parallel for collapse(2)
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    datResize[ii][t][k] = datDS[t][k];  //zScoreMap scaling
                }//for(k)
            }//for(t)

//            datResize[ii].resize(T);
//            for (int t = 0; t < T; ++t) {
//                datResize[ii][t].resize(L);
//                for (int k = 0; k < L; ++k) {
//                    datResize[ii][t][k] = datDS[t][k];
//                }//for(k)
//            }//for(t)


        }//for(ii)
//        cout<<datResize.size()<<endl;
//        cout<<datResize[0].size()<<endl;
//        for (int i = 0; i < 5; ++i) {
//            for (int j = 0; j < 5; ++j) {
//                cout<<datResize[0][0][0].at<float>(i,j)<<" ";
//            }
//            cout<<endl;
//        }
        return datResize;
    }//normalizeAndResize()
    //                !!!!!!!!!!!!!!!!LAST CORRECT!!!!!!!!!!!!!!!!!!!!!


    void seedDetect2_DS_accelerate(vector<vector<cv::Mat>> dF, const vector<vector<cv::Mat>>& dataOrg,
                                   const vector<vector<int>>& arLst){

        vector<float> Thrs;
        for (float i = opts.maxdF1; i >=opts.thrARScl ; i=i-opts.step) {
            Thrs.emplace_back(i);
        }
        vector<float> scaleRatios = {2,4,8};
        opts.scaleRatios = scaleRatios;
        int H = dF[0][0].rows;
        int W = dF[0][0].cols;
        int L = dF[0].size();
        int T = dF.size();

        //assign saturation part can always be selected when checking seed
        if (opts.maxValueDat1 == pow(2, opts.BitDepth) - 1){
//            #pragma omp parallel for collapse(2)
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    cv::Mat mask = (dataOrg[t][k]==1);
                    dF[t][k].setTo(INFINITY,mask);
                }
            }
        }//if

        vector<int> regSz(arLst.size());
        vector<vector<cv::Mat>> activeMap(T,vector<cv::Mat>(L));
//        #pragma omp parallel for collapse(2)
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                activeMap[t][k] = cv::Mat::zeros(H,W,CV_16U);
            }
        }
        //active map
        for (int ii = 0; ii < arLst.size(); ++ii) { // ii -- different regions
            unordered_set<float> unique_indices;
            for (int ii_ite = 0; ii_ite < arLst[ii].size(); ++ii_ite) { //ii_ite -- different points in a region
                int i = ind2sub(arLst[ii][ii_ite]-1,H,W,L).i; // WARNING: using ind from matlab, do not -1
                int j = ind2sub(arLst[ii][ii_ite]-1,H,W,L).j;
                int k = ind2sub(arLst[ii][ii_ite]-1,H,W,L).k;
                int t = ind2sub(arLst[ii][ii_ite]-1,H,W,L).t;
                activeMap[t][k].at<ushort>(i,j) = ii+1;
                
                int ind = sub2ind(i,j,k,H,W);
                unique_indices.insert(ind);
            }//for(ii_ite)
            regSz[ii] = unique_indices.size();
//            cout<<regSz[ii]<<" ";
        }//for(ii)

        //down sampled data
        vector<vector<vector<cv::Mat>>> validMaps(scaleRatios.size(),vector<vector<cv::Mat>>(T,vector<cv::Mat>(L))); //treat as boolean
        /*
         * datResize is flattened in MATLAB, remain 4D matrix in C++
         */
        vector<vector<vector<cv::Mat>>> datResize = normalizeAndResize(dataOrg); //normalized data to do significance test
        vector<vector<vector<cv::Mat>>> dFResize(scaleRatios.size(),vector<vector<cv::Mat>>(T,vector<cv::Mat>(L))); //down sampled data to do selection
        vector<float> H0s(scaleRatios.size(), 0);
        vector<float> W0s(scaleRatios.size(), 0);
        for (int j = 0; j < scaleRatios.size(); ++j) {
//            datResize[j] = reshape
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
//                    cv::Mat temp_dF;
//                    cv::Mat temp_validMaps;
//                    cv::resize(dF[t][k], temp_dF, cv::Size(), 1 / scaleRatios[j], 1 / scaleRatios[j], cv::INTER_AREA);
//                    cv::resize(activeMap[t][k], temp_validMaps, cv::Size(), 1 / scaleRatios[j], 1 / scaleRatios[j], cv::INTER_AREA);
//                    temp_dF = myResize(dF[t][k], scaleRatios[j], scaleRatios[j]);
//                    temp_validMaps = myResize(activeMap[t][k],scaleRatios[j],scaleRatios[j]);
                    dFResize[j][t][k] = myResize(dF[t][k], scaleRatios[j], scaleRatios[j]);
                    validMaps[j][t][k] = myResize(activeMap[t][k],scaleRatios[j],scaleRatios[j]);
                }//for(k)
            }//for(t)
            H0s[j] = ceil(H/scaleRatios[j]);
            W0s[j] = ceil(W/scaleRatios[j]);
        }//for(j)

        //seed map
        vector<vector<cv::Mat>> zscoreMap(T,vector<cv::Mat>(L));
//        #pragma omp parallel for collapse(2)
        for (int t = 0; t < dF.size(); ++t) {
            for (int k = 0; k < dF[0].size(); ++k) {
                zscoreMap[t][k] = cv::Mat::zeros(dF[0][0].rows,dF[0][0].cols, CV_32F);
            }
        }//for(t)

        for (int ii = 0; ii < Thrs.size(); ++ii) { // threshold --k
            float curThr = Thrs[ii];

            for (int ii_ds = 0; ii_ds < scaleRatios.size(); ++ii_ds) { // down sample rate --j
                float H0 = H0s[ii_ds];
                float W0 = W0s[ii_ds];
                float scaleRatio = scaleRatios[ii_ds];
                vector<int> tmp;

                for (int ind = 0; ind < scaleRatio; ++ind) {
                    for (int ind_i = 0; ind_i < scaleRatio; ++ind_i) {
                        tmp.emplace_back(ind+1);
                    }
                }//for(ind)

                vector<vector<cv::Mat>> selectMap(dFResize[ii_ds].size(),vector<cv::Mat>(dFResize[ii_ds][0].size()));
//                #pragma omp parallel for collapse(2)
                for (int t = 0; t < dFResize[ii_ds].size(); ++t) {
                    for (int k = 0; k < dFResize[ii_ds][0].size(); ++k) {
                        selectMap[t][k] = cv::Mat::zeros(dFResize[ii_ds][0][0].rows, dFResize[ii_ds][0][0].cols, CV_8U);
                        cv::Mat mask = ((dFResize[ii_ds][t][k]>curThr) & (validMaps[ii_ds][t][k]!=0));
                        selectMap[t][k].setTo(1,mask);
//                        cv::Mat nonZero;
//                        cv::findNonZero(selectMap[t][k],nonZero);
//                        for (int i = 0; i < nonZero.total(); ++i) {
//                            cout<<sub2ind(nonZero.at<cv::Point>(i).x,nonZero.at<cv::Point>(i).y,k,t,
//                                          selectMap[0][0].rows,selectMap[0][0].cols,dFResize[ii_ds][0].size())<<" ";
//                        }
//                        if (cv::countNonZero(selectMap[t][k])>0){
//                            cout<<"error"<< endl;
//                        }

//                        for (int i = 0; i < dFResize[0][0][0].rows; ++i) {
//                            for (int j = 0; j < dFResize[0][0][0].cols; ++j) {
//                                if ((dFResize[ii_ds][t][k].at<float>(i,j) > curThr) && (validMaps[ii][t][k].at<int>(i,j) !=0)){
//                                    selectMap[t][k].at<unsigned char>(i,j) = 1;
//                                }//if
//                            }//for(j)
//                        }//for(i)

                    }//for(k)
                }//for(t)
//                writeDataToMatFile(selectMap, "C:/Users/Kevin Qiao/Desktop/AQuA_data/test/sel.mat");
                vector<vector<int>> curRegions = bw2Reg(selectMap);
                // rough filter -- for acceleration
//                for (int ii_cur = 0; ii_cur < curRegions.size(); ++ii_cur) {
//                    if (curRegions[ii_cur].size() <= (opts.minSize / pow(scaleRatio,2) * opts.minDur / 3)){
//                        curRegions.erase(curRegions.begin() + ii_cur);
//                        --ii_cur;
//                    }
//                }
                curRegions.erase(remove_if(curRegions.begin(), curRegions.end(),
                    [&](const auto& region) {
                        return region.size() <= (opts.minSize / pow(scaleRatio, 2) * opts.minDur / 3);
                    }), curRegions.end()); //move all the wanted elements to the front, and iterator points at next pos of the element wanted, then delete


                for (int ii_cur = 0; ii_cur < curRegions.size(); ++ii_cur) { //--i
                    unordered_set<int> ihw_temp;
                    vector<int> pix;
                    vector<int> ih;
                    vector<int> iw;
                    vector<int> il;
                    vector<int> it;
                    for (int ii_ite = 0; ii_ite < curRegions[ii_cur].size(); ++ii_ite) {
                        pix.emplace_back(curRegions[ii_cur][ii_ite]);
                    }
                    sort(pix.begin(), pix.end());

                    for (int ii_ite = 0; ii_ite < pix.size(); ++ii_ite) {
                        Point_struct pix_temp = ind2sub(pix[ii_ite],H0,W0,L);
                        ih.emplace_back(pix_temp.i);
                        iw.emplace_back(pix_temp.j);
                        il.emplace_back(pix_temp.k);
                        it.emplace_back(pix_temp.t);
                        ihw_temp.insert(sub2ind(ih[ii_ite],iw[ii_ite],il[ii_ite],H0,W0));
                    }
                    vector<int> ihw(ihw_temp.begin(), ihw_temp.end());
                    sort(ihw.begin(), ihw.end());
                    int dur = *max_element(it.begin(),it.end()) - *min_element(it.begin(), it.end()) + 1;
                    int arLabel_hs = ih[0] * scaleRatio;
                    int arLabel_he = min(static_cast<float>(H),ih[0]*scaleRatio + 1);
                    int arLabel_ws = iw[0] * scaleRatio;
                    int arLabel_we = min(static_cast<float>(W),iw[0]*scaleRatio + 1);
                    vector<int> arLabel;
                    for (int i = arLabel_hs; i <= arLabel_he; ++i) {
                        for (int j = arLabel_ws; j <= arLabel_we; ++j) {
                            if (activeMap[it[0]][il[0]].at<ushort>(i, j) != 0) {
                                arLabel.emplace_back(static_cast<int>(activeMap[it[0]][il[0]].at<ushort>(i, j)));
                            }
                        }
                    }
                    int arLabel_val = *min_element(arLabel.begin(), arLabel.end());

                    //filter according to size and duration, also check seed detected or not
                    if (dur < opts.minDur || ihw.size()< max(static_cast<float>(opts.minSize), regSz[arLabel_val-1] * opts.seedSzRatio) /
                                                         pow(scaleRatio,2)){
                        continue;
                    }

                    // convert back
                    cv::Mat ihOrg = cv::Mat(ih.size(), scaleRatio*scaleRatio, CV_16U);
                    cv::Mat iwOrg = cv::Mat(ih.size(), scaleRatio*scaleRatio, CV_16U);
                    cv::Mat ilOrg = cv::Mat(ih.size(), scaleRatio*scaleRatio, CV_16U);
                    cv::Mat itOrg = cv::Mat(ih.size(), scaleRatio*scaleRatio, CV_16U);
                    cv::Mat select = cv::Mat::zeros(ih.size(), scaleRatio*scaleRatio, CV_8U);

                    for (int col = 0, add = 1; col < pow(scaleRatio,2); ++col,++add) {
                        if (add == scaleRatio+1){
                            add = 1;
                        }
//                        #pragma omp parallel for
                        for (int row = 0; row < ih.size(); ++row) {
                            ihOrg.at<ushort>(row,col) = ih[row]*scaleRatio - 1 + add;  // matlab: (ih[]-1)*scaleRatio
                            iwOrg.at<ushort>(row,col) = iw[row]*scaleRatio - 1 + tmp[col];
                            ilOrg.at<ushort>(row,col) = il[row];
                            itOrg.at<ushort>(row,col) = it[row];
                        }
                    }//for(col)

//                    cv::Mat mask = ((ihOrg < H) & (iwOrg < W));
                    select.setTo(1,((ihOrg < H) & (iwOrg < W)));
//                    for (int i = 0; i < ih.size(); ++i) {
//                        for (int j = 0; j < scaleRatio * scaleRatio; ++j) {
//                            if ((ihOrg.at<float>(i,j) <= H) && (iwOrg.at<float>(i,j) <= W)){
//                                select.at<uint>(i,j) = 1;
//                            }
//                        }
//                    }

                    vector<int> pixOrg;
                    for (int j = 0; j < scaleRatio * scaleRatio; ++j) {
                        for (int i = 0; i < ih.size(); ++i) {
                            if (select.at<uchar>(i,j) == 1) {
                                pixOrg.emplace_back(sub2ind(ihOrg.at<ushort>(i,j),iwOrg.at<ushort>(i,j),
                                                            ilOrg.at<ushort>(i,j),itOrg.at<ushort>(i,j),H,W,L));
                            }
                        }
                    }

                    bool notEmpty = false;
                    for (int num = 0; num < pixOrg.size(); ++num) {
                        Point_struct pix_temp = ind2sub(pixOrg[num],H,W,L);
                        if (zscoreMap[pix_temp.t][pix_temp.k].at<float>(pix_temp.i,pix_temp.j)>0) {
                            notEmpty = true;
                            break;
                        }
                    }

                    if ((pixOrg.size()< max(static_cast<float>(opts.minSize), regSz[arLabel_val-1]*opts.seedSzRatio)) || notEmpty){
                        continue;
                    }

                    //calculate significance
                    float t_scl = max(1.0, round(static_cast<double>(dur)/opts.TPatch));
                    getSeedScore_DS4(pix,datResize[ii_ds],H0,W0,L,T,t_scl);




                }//for(ii_cur) --curRegions --i
            }//for(ii_ds) --scaleRatios --j
        }//for(ii) --Thrs --k
        




    }//seedDetect2_DS_accelerate


    void seDetection(const vector<vector<cv::Mat>>& dF, const vector<vector<cv::Mat>>& dataOrg,
                     const vector<vector<int>>& arLst){
        int H = dF[0][0].rows;
        int W = dF[0][0].cols;
        int L = dF[0].size();
        int T = dF.size();
        //setting
        //may be according to active region size, set different scale
        for (int i = opts.maxSpaScale; i >= opts.minSpaScale; --i) {
            opts.scaleRatios.emplace_back(i);
        }
        cout<< "--------start seed detecting--------"<<endl;
        /*
         * test split by xuelong's methods from active region
         * top-down seed detection
         * is there any advantage that when checking seeds, calculate the score
         * for each pixel then combine them together? -- Since each individual
         * pixel in seed has different time windows. Consider the potential
         * propagation, we cannot use the average curve for seed detection.
         */
        seedDetect2_DS_accelerate(dF,dataOrg,arLst);




    }//seDetection()


    void phaseRun(){
        cout<< "--------start temporal segmentation--------"<<endl;
        if (opts.needTemp){
            opts.step = 0.5; // sigma
            opts.sigThr = 3.5; // Z score of seed significance
            opts.maxDelay = 0.6; // allowed maximum dissimilarity in merging
            opts.seedSzRatio = 0.01; // seed size relative to active region
            opts.needRefine = false; // peaks are temporally adjacent, need regine
            opts.needGrow = false; // grow active regions according to signal pattern
            opts.maxdF1 = 15.4560;

//            vector<vector<cv::Mat>> dF1 = opts.dF1;
//            vector<vector<cv::Mat>> dataOrg1 = opts.data1_org;
//            vector<vector<Point_struct>> arLst1 = opts.arLst1;
            vector<vector<cv::Mat>> dataOrg1 = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/phaseRun.mat", "datOrg1");
            vector<vector<cv::Mat>> dF1 = AQuA::load4D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/phaseRun.mat", "dF1");
            vector<vector<int>> arLst1 = AQuA::loadCell_int("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/phaseRun.mat", "arLst1");
            AQuA::opts.tempVarOrg1 = AQuA::load3D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/tempVar.mat", "tempVar");
            AQuA::opts.correctPars1 = AQuA::load3D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/tempVar.mat", "correctPars");
            seDetection(dF1,dataOrg1,arLst1);



        }//if

    }//phaseRun


}//namespace
