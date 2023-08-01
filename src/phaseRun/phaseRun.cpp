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

/*
    Score_struct getSeedScore_DS4(vector<Point_struct> pix, vector<vector<cv::Mat>> datVec, int H0, int W0, int t_scl){
        Score_struct result;
        //down sample
        for (int i = 0; i < pix.size(); ++i) {
            pix[i].t = pix[i].t/t_scl;
        }
        int T0 = floor(T/t_scl);
        cv::Mat select = cv::Mat::zeros(pix.size(),1,CV_8U);
        for (int i = 0; i < pix.size(); ++i) {
            if (pix[i].t <= T0){
                select.at<uint>(i,0) = 1;
            }
        }
        unordered_set<int> fgPix;
        unordered_map<int,int> a_counts;
        for (int i = 0; i < pix.size(); ++i) {
            if (select.at<uint>(i,0) == 1){
                int ind = sub2ind(pix[i].i,pix[i].j,pix[i].k,pix[i].t,H0,W0,L);
                a_counts[ind]++;
            }
        }
        for (const auto& pair: a_counts) {
            if (pair.second > (t_scl/2)){
                fgPix.insert(pair.first);
            }
        }

        if (fgPix.size() ==0){
            result.z_score1 = 0;
            result.z_score2 = 0;
            result.t_score1 = 0;
            result.t_score2 = 0;
            return result;
        }

        //find neighbor
        vector<int> ih;
        vector<int> iw;
        vector<int> il;
        vector<int> it;
        for (auto const& fg:fgPix) {
            int h = ind2sub(fg,H0,W0,L).i;
            int w = ind2sub(fg,H0,W0,L).j;
            int k = ind2sub(fg,H0,W0,L).k;
            int t = ind2sub(fg,H0,W0,L).t;
            ih.emplace_back(h);
            iw.emplace_back(w);
            il.emplace_back(k);
            it.emplace_back(t);
        }
        vector<int> ihwOrg;
        unordered_set<int> ihw;
        for (int i = 0; i < ih.size(); ++i) {
            ihwOrg.emplace_back(sub2ind(ih[i],iw[i],il[i],H0,W0));
            ihw.insert(ihwOrg[i]);
        }
        vector<int> sz;
        vector<vector<float>> fgAll(ihw.size());
        vector<vector<float>> bgL(ihw.size());
        vector<vector<float>> bgR(ihw.size());
        vector<vector<float>> nanVec(ihw.size());
        int cnt = 0;
        int cnt2 = 0;
        int count = 0;
        vector<float> noise;
        vector<int> degreeOfFreedoms;
        for (const auto& i:ihw) {
            vector<int> curIt;
            for (int ii = 0; ii < ihwOrg.size(); ++ii) {
                if (i == ihwOrg[ii]){
                    curIt.emplace_back(it[ii]);
                }
            }//for(ii)
            //get normalized down sampled curve
            cv::Mat curve0 = cv::Mat(T,1,CV_32F);
            int h_ihw = ind2sub(i,H0,W0).i;
            int w_ihw = ind2sub(i,H0,W0).j;
            int k_ihw = ind2sub(i,H0,W0).k;
            for (int t = 0; t < T; ++t) {
                curve0.at<float>(t,0) = datVec[t][k_ihw].at<float>(h_ihw,w_ihw);
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
            for (const auto& ite:curIt) {
                fgAll[count].emplace_back(curve.at<float>(ite,0));
            }
            for (const auto& ite:t_Left) {
                bgL[count].emplace_back(curve.at<float>(ite,0));
            }
            for (const auto& ite:t_Right) {
                bgR[count].emplace_back(curve.at<float>(ite,0));
            }
            vector<int> diff;
            vector<int> t0_t1;
            for (int ii = t0; ii <= t1; ++ii) {
                t0_t1.emplace_back(ii);
            }
            set_difference(t0_t1.begin(),t0_t1.end(),curIt.begin(),curIt.end(),
                                inserter(diff,diff.begin()));
            for (const auto& ite:diff) {
                nanVec[count].emplace_back(curve.at<float>(ite,0));
            }

            //jump one point
            vector<float> difL;
            vector<float> difR;
            if (!t_Left.empty()){
                for (auto& ite: t_Left) {
                    ite -= 1;
                }
                t0 = max(1,*min_element(t_Left.begin(),t_Left.end())*t_scl);
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
                t1 = min(*max_element(t_Right.begin(),t_Right.end())*t_scl, T-2);
                for (int tt = t0;tt <= t1;++tt) {
                    difR.emplace_back(pow((curve0.at<float>(tt,0)-curve0.at<float>(tt+1,0)),2));
                }
                for (int tt = t0; tt <= t1+1; ++tt) {
                    if(curve0.at<float>(tt,0) == 0){
                        ++cnt;
                    }
                }
                cnt2 = cnt2 + t1 -t0 + 2;
            }//if(!t_left.empty())
            noise.insert(noise.end(),difL.begin(),difL.end());
            noise.insert(noise.end(),difR.begin(),difR.end());
            degreeOfFreedoms.emplace_back(difL.size());
            degreeOfFreedoms.emplace_back(difR.size());
            count++;
        }//for(ihw)




    }//getSeedScore_DS4
    */


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
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
//                    cv::resize(dataOrg[t][k],datDS[t][k],cv::Size(),1/scaleRatio,1/scaleRatio, cv::INTER_AREA);
                    datDS[t][k] = myResize(dataOrg[t][k],scaleRatio,scaleRatio);
                }//for(k)
            }//for(t)

            //consider the possible noise correlation, need to re-estimate noise
            vector<cv::Mat> curVarMap(L);
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
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    cv::Mat mask = (dataOrg[t][k]==1);
                    dF[t][k].setTo(INFINITY,mask);
                }
            }
        }//if

        vector<int> regSz(arLst.size());
        vector<vector<cv::Mat>> activeMap(T,vector<cv::Mat>(L));
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
        vector<vector<vector<cv::Mat>>> datResize = normalizeAndResize(dataOrg); //normalized data to do significance test
        vector<vector<vector<cv::Mat>>> dFResize(scaleRatios.size(),vector<vector<cv::Mat>>(T,vector<cv::Mat>(L))); //down sampled data to do selection
        vector<float> H0s(scaleRatios.size(), 0);
        vector<float> W0s(scaleRatios.size(), 0);
        for (int ii = 0; ii < scaleRatios.size(); ++ii) {
//            datResize[ii] = reshape
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
//                    cv::Mat temp_dF;
//                    cv::Mat temp_validMaps;
//                    cv::resize(dF[t][k], temp_dF, cv::Size(), 1 / scaleRatios[ii], 1 / scaleRatios[ii], cv::INTER_AREA);
//                    cv::resize(activeMap[t][k], temp_validMaps, cv::Size(), 1 / scaleRatios[ii], 1 / scaleRatios[ii], cv::INTER_AREA);
//                    temp_dF = myResize(dF[t][k], scaleRatios[ii], scaleRatios[ii]);
//                    temp_validMaps = myResize(activeMap[t][k],scaleRatios[ii],scaleRatios[ii]);
                    dFResize[ii][t][k] = myResize(dF[t][k], scaleRatios[ii], scaleRatios[ii]);
                    validMaps[ii][t][k] = myResize(activeMap[t][k],scaleRatios[ii],scaleRatios[ii]);
                }//for(k)
            }//for(t)
            H0s[ii] = ceil(H/scaleRatios[ii]);
            W0s[ii] = ceil(W/scaleRatios[ii]);
        }//for(ii)

        //seed map
        vector<vector<cv::Mat>> zscoreMap(T,vector<cv::Mat>(L));
        for (int t = 0; t < dF.size(); ++t) {
            for (int k = 0; k < dF[0].size(); ++k) {
                zscoreMap[t][k] = cv::Mat::zeros(dF[0][0].rows,dF[0][0].cols, CV_32F);
            }
        }//for(t)

        for (int ii = 0; ii < Thrs.size(); ++ii) { // threshold
            float curThr = Thrs[ii];

            for (int ii_ds = 0; ii_ds < scaleRatios.size(); ++ii_ds) { // down sample rate
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
                for (int t = 0; t < dFResize[ii_ds].size(); ++t) {

                    for (int k = 0; k < dFResize[ii_ds][0].size(); ++k) {
                        selectMap[t][k] = cv::Mat::zeros(dFResize[ii_ds][0][0].rows, dFResize[ii_ds][0][0].cols, CV_8U);
                        cv::Mat mask = ((dFResize[ii_ds][t][k]>curThr) & (validMaps[ii_ds][t][k]!=0));
                        selectMap[t][k].setTo(1,mask);
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
                for (int ii_cur = 0; ii_cur < curRegions.size(); ++ii_cur) {
                    if (curRegions[ii_cur].size() <= (opts.minSize / pow(scaleRatio,2) * opts.minDur / 3)){
                        curRegions.erase(curRegions.begin() + ii_cur);
                    }
                }

                for (int ii_cur = 0; ii_cur < curRegions.size(); ++ii_cur) {
                    unordered_set<int> ihw;
                    vector<int> ih;
                    vector<int> iw;
                    vector<int> il;
                    vector<int> it;
                    for (int ii_ite = 0; ii_ite < curRegions[ii_cur].size(); ++ii_ite) {
                        int pix = curRegions[ii_cur][ii_ite];

                        ih.emplace_back(ind2sub(pix,H0,W0,L).i);
                        iw.emplace_back(ind2sub(pix,H0,W0,L).j);
                        il.emplace_back(ind2sub(pix,H0,W0,L).k);
                        it.emplace_back(ind2sub(pix,H0,W0,L).t);
                        ihw.insert(sub2ind(ih[ii_ite],iw[ii_ite],il[ii_ite],H0,W0));
                    }

                    int dur = *max_element(it.begin(),it.end()) - *min_element(it.begin(), it.end());
                    int arLabel_hs = ih[0] * scaleRatio;
                    int arLabel_he = min(static_cast<float>(H),ih[0]*scaleRatio + 1);
                    int arLabel_ws = iw[0] * scaleRatio;
                    int arLabel_we = min(static_cast<float>(W),iw[0]*scaleRatio + 1);
                    vector<int> arLabel;
                    for (int i = arLabel_hs; i <= arLabel_he; ++i) {
                        for (int j = arLabel_ws; j < arLabel_we; ++j) {
                            if (activeMap[it[0]][il[0]].at<ushort>(i, j) != 0) {
                                arLabel.emplace_back(static_cast<int>(activeMap[it[0]][il[0]].at<ushort>(i, j)));
                            }
                        }
                    }
                    int arLabel_val = *min_element(arLabel.begin(), arLabel.end());

                    //filter according to size and duration, also check seed detected or not
                    if (dur < opts.minDur || ihw.size()< max(static_cast<float>(opts.minSize), regSz[arLabel_val]*opts.seedSzRatio) /
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

                        for (int row = 0; row < ih.size(); ++row) {
                            ihOrg.at<ushort>(row,col) = ih[row]*scaleRatio + add;
                            iwOrg.at<ushort>(row,col) = iw[row]*scaleRatio + add;
                            ilOrg.at<ushort>(row,col) = il[row];
                            itOrg.at<ushort>(row,col) = it[row];
                        }
                    }//for(col)

                    for (int i = 0; i < ih.size(); ++i) {
                        for (int j = 0; j < scaleRatio * scaleRatio; ++j) {
                            if ((ihOrg.at<float>(i,j) <= H) && (iwOrg.at<float>(i,j) <= W)){
                                select.at<uint>(i,j) = 1;
                            }
                        }
                    }
                    vector<Point_struct> pixOrg;
                    for (int i = 0; i < ih.size(); ++i) {
                        for (int j = 0; j < scaleRatio * scaleRatio; ++j) {
                            if (select.at<uint>(i,j) == 1) {
                                Point_struct temp;
                                temp.i = ihOrg.at<float>(i,j);
                                temp.j = iwOrg.at<float>(i,j);
                                temp.k = ilOrg.at<float>(i,j);
                                temp.t = itOrg.at<float>(i,j);
                                pixOrg.emplace_back(temp);
                            }
                        }
                    }
                    bool find_flag = false;
                    for (int num = 0; num < pixOrg.size(); ++num) {
                        int i = pixOrg[num].i;
                        int j = pixOrg[num].j;
                        int k = pixOrg[num].k;
                        int t = pixOrg[num].t;
                        if (zscoreMap[t][k].at<float>(i,j)>0){
                            find_flag = true;
                            break;
                        }
                    }
                    if ((pixOrg.size()< max(static_cast<float>(opts.minSize), regSz[arLabel_val]*opts.seedSzRatio)) && find_flag){
                        continue;
                    }

                    //calculate significance
                    float t_scl =  max(static_cast<double>(1), round(dur/opts.TPatch));



                }//for(ii_cur)
            }//for(ii_ds) --scaleRatios.size()
        }//for(ii) --Thrs.size()
        




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
            vector<vector<int>> arLst1 = AQuA::loadCell("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/phaseRun.mat", "arLst1");
            AQuA::opts.tempVarOrg1 = AQuA::load3D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/tempVar.mat", "tempVar");
            AQuA::opts.correctPars1 = AQuA::load3D("C:/Users/Kevin Qiao/Desktop/AQuA_data/test/tempVar.mat", "correctPars");
            seDetection(dF1,dataOrg1,arLst1);



        }//if

    }//phaseRun


}//namespace
