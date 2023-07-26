//
// Created by Kevin Qiao on 7/18/2023.
//

#include "phaseRun.h"

namespace AQuA{

    cv::Mat myResize(const cv::Mat& src, float scaleRatio_x, float scaleRatio_y){
        float h = std::ceil(src.rows/scaleRatio_y); //height of dst
        float w = std::ceil(src.cols/scaleRatio_x); // width of dst
        int unit_h = std::ceil(src.rows / h); //number of units in row
        int unit_w = std::ceil(src.cols / w); // number of units in col
        cv::Mat dst = cv::Mat(h, w, CV_32F);
        for (int i_dst = 0; i_dst < h; ++i_dst) {
            for (int j_dst = 0; j_dst < w; ++j_dst) {
                float sum = 0;
                float count = 0;
                for (int i_src = i_dst*unit_h; (i_src<(i_dst+1)*unit_h) && (i_src<src.rows); ++i_src) {
                    for (int j_src = j_dst*unit_w; (j_src<(j_dst+1)*unit_w) && (j_src<src.cols); ++j_src) {
                        sum += src.at<float>(i_src,j_src);
                        ++count;
//                        std::cout<<"i_src: "<<i_src<<"  j_src: "<<j_src<<" ";
                    }//for(j_src)
//                    std::cout<<std::endl;
                }//for(i_src)
                dst.at<float>(i_dst, j_dst) = sum/count;
            }
        }
        return dst;
    }//myResize()


    Score_struct getSeedScore_DS4(std::vector<Point_struct> pix, std::vector<std::vector<cv::Mat>> datVec, int H0, int W0, int t_scl){
        Score_struct result;
        //down sample
        for (int i = 0; i < pix.size(); ++i) {
            pix[i].t = pix[i].t/t_scl;
        }
        int T0 = std::floor(T/t_scl);
        cv::Mat select = cv::Mat::zeros(pix.size(),1,CV_8U);
        for (int i = 0; i < pix.size(); ++i) {
            if (pix[i].t <= T0){
                select.at<uint>(i,0) = 1;
            }
        }
        std::unordered_set<int> fgPix;
        std::unordered_map<int,int> a_counts;
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
        std::vector<int> ih;
        std::vector<int> iw;
        std::vector<int> il;
        std::vector<int> it;
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
        std::vector<int> ihwOrg;
        std::unordered_set<int> ihw;
        for (int i = 0; i < ih.size(); ++i) {
            ihwOrg.emplace_back(sub2ind(ih[i],iw[i],il[i],H0,W0));
            ihw.insert(ihwOrg[i]);
        }
        std::vector<int> sz;
        std::vector<std::vector<float>> fgAll(ihw.size());
        std::vector<std::vector<float>> bgL(ihw.size());
        std::vector<std::vector<float>> bgR(ihw.size());
        std::vector<std::vector<float>> nanVec(ihw.size());
        int cnt = 0;
        int cnt2 = 0;
        int count = 0;
        for (const auto& i:ihw) {
            std::vector<int> curIt;
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
            cv::multiply(curve,std::sqrt(t_scl),curve);
            int dur = curIt.size();
            sz.emplace_back(dur);
            int t0 = *std::min_element(curIt.begin(),curIt.end());
            int t1 = *std::max_element(curIt.begin(),curIt.end());
            std::vector<int> t_Left;
            std::vector<int> t_Right;
            for (int tt = std::max(0,t0-dur); tt <= t0 -1; ++tt) {
                t_Left.emplace_back(tt);
            }
            for (int tt = t1+1; tt <= std::min(T0, t1 + dur); ++tt) {
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
            std::vector<int> diff;
            std::vector<int> t0_t1;
            for (int ii = t0; ii <= t1; ++ii) {
                t0_t1.emplace_back(ii);
            }
            std::set_difference(t0_t1.begin(),t0_t1.end(),curIt.begin(),curIt.end(),
                                std::inserter(diff,diff.begin()));
            for (const auto& ite:diff) {
                nanVec[count].emplace_back(curve.at<float>(ite,0));
            }

            //jump one point
            std::vector<float> difL;
            if (!t_Left.empty()){
                for (auto& ite: t_Left) {
                    ite -= 1;
                }
                t0 = std::max(2,*std::max_element(t_Left.begin(),t_Left.end())*t_scl);
                t1 = *std::max_element(t_Left.begin(),t_Left.end())*t_scl;
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




            count++;
        }//for(ihw)




    }//getSeedScore_DS4


    std::vector<std::vector<std::vector<cv::Mat>>> normalizeAndResize(const std::vector<std::vector<cv::Mat>>& dataOrg){
        std::vector<float> scaleRatios = opts.scaleRatios;
        std::vector<std::vector<cv::Mat>> datDS(T);
        std::vector<std::vector<std::vector<cv::Mat>>> datResize(scaleRatios.size());

        //rescl dForg
        for (int ii = 0; ii < scaleRatios.size(); ++ii) {
            int scaleRatio = scaleRatios[ii];
            for (int t = 0; t < T; ++t) {
                datDS[t].resize(L);
                for (int k = 0; k < L; ++k) {
//                    cv::resize(dataOrg[t][k],datDS[t][k],cv::Size(),1/scaleRatio,1/scaleRatio, cv::INTER_AREA);
                    datDS[t][k] = myResize(dataOrg[t][k],scaleRatio,scaleRatio);
                }//for(k)
            }//for(t)
            //consider the possible noise correlation, need to re-estimate noise
            std::vector<cv::Mat> curVarMap(L);
            for (int k = 0; k < L; ++k) {
                curVarMap.emplace_back(cv::Mat(datDS[0][0].rows,datDS[0][0].cols,CV_32F));
                for (int i = 0; i < datDS[0][0].rows; ++i) {
                    for (int j = 0; j < datDS[0][0].cols; ++j) {
                        double sum = 0;
                        for (int t = 1; t < T; ++t) {
                            sum += pow((datDS[t][k].at<float>(i,j) - datDS[t-1][k].at<float>(i,j)), 2);
                        }
                        curVarMap[k].at<float>(i,j) = static_cast<float>(sum/T/2);
                    }
                }
            }//for(k)
            //correct noise due to truncation
            std::vector<cv::Mat> var1(L);
            std::vector<cv::Mat> var2(L);
            std::vector<cv::Mat> curStdMap(L);
            for (int k = 0; k < L; ++k) {
                cv::Mat temp1;
                cv::multiply(opts.tempVarOrg1[k], 2, temp1);
                cv::divide(temp1, opts.correctPars1[k], temp1);
//                cv::resize(opts.tempVarOrg1[k], var1[k], cv::Size(), 1 / scaleRatio, 1 / scaleRatio, cv::INTER_AREA);
//                cv::resize(temp1, var2[k], cv::Size(), 1 / scaleRatio, 1 / scaleRatio, cv::INTER_AREA);
                var1[k] = myResize(opts.tempVarOrg1[k],scaleRatio,scaleRatio);
                var2[k] = myResize(temp1, scaleRatio,scaleRatio);


                cv::Mat temp2;
                cv::multiply(curVarMap[k], var2[k], temp2);
                cv::divide(temp2, var1[k], temp2);
                cv::sqrt(temp2, curStdMap[k]);
                cv::divide(datDS[k], curStdMap[k], datDS[k]); // since it is only used in seed detection,
                // and seed detection only check single pixel's score. Not need to fitting again
            }//for(k)

            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    datResize[ii][t].emplace_back(datDS[t][k]);
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

        return datResize;
    }//normalizeAndResize()


    void seedDetect2_DS_accelerate(std::vector<std::vector<cv::Mat>> dF, const std::vector<std::vector<cv::Mat>>& dataOrg,
                                   std::vector<std::vector<Point_struct>>& arLst){

        std::vector<float> Thrs;
        for (int i = opts.maxdF1; i >=opts.thrARScl ; i=i-opts.step) {
            Thrs.emplace_back(i);
        }
        std::vector<float> scaleRatios = {2,4,8};
        opts.scaleRatios = scaleRatios;

        //assign saturation part can always be selected when checking seed
        if (opts.maxValueDat1 == pow(2, opts.BitDepth) - 1){
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    for (int i = 0; i < H; ++i) {
                        for (int j = 0; j < W; ++j) {
                            if (dataOrg[t][k].at<float>(i,j) == 1){
                                dF[t][k].at<float>(i,j) = INFINITY;
                            }
                        }
                    }
                }
            }
        }//if
        std::vector<float> regSz(arLst.size());
        std::vector<std::vector<cv::Mat>> activeMap(T);
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                activeMap[t].emplace_back(cv::Mat::zeros(H,W,CV_32S));
            }
        }
        //active map
        for (int ii = 0; ii < arLst.size(); ++ii) { // ii -- different regions
            std::unordered_set<float> unique_indices;
            for (int ii_ite = 0; ii_ite < arLst[ii].size(); ++ii_ite) { //ii_ite -- different points in a region
                int i = arLst[ii][ii_ite].i;
                int j = arLst[ii][ii_ite].j;
                int k = arLst[ii][ii_ite].k;
                int t = arLst[ii][ii_ite].t;
                activeMap[t][k].at<int>(i,j) = ii+1;
                int ind = sub2ind(i,j,k,H,W);
                unique_indices.insert(ind);
            }//for(ii_ite)
            regSz[ii] = unique_indices.size();
        }//for(ii)

        //down sampled data
        std::vector<std::vector<std::vector<cv::Mat>>> validMaps(scaleRatios.size()); //treat as boolean
        std::vector<std::vector<std::vector<cv::Mat>>> datResize = normalizeAndResize(dataOrg); //normalized data to do significance test
        std::vector<std::vector<std::vector<cv::Mat>>> dFResize(scaleRatios.size()); //down sampled data to do selection
        std::vector<float> H0s(scaleRatios.size(), 0);
        std::vector<float> W0s(scaleRatios.size(), 0);
        for (int ii = 0; ii < scaleRatios.size(); ++ii) {
//            datResize[ii] = reshape
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    cv::Mat temp_dF;
                    cv::Mat temp_validMaps;
//                    cv::resize(dF[t][k], temp_dF, cv::Size(), 1 / scaleRatios[ii], 1 / scaleRatios[ii], cv::INTER_AREA);
//                    cv::resize(activeMap[t][k], temp_validMaps, cv::Size(), 1 / scaleRatios[ii], 1 / scaleRatios[ii], cv::INTER_AREA);
                    temp_dF = myResize(dF[t][k], scaleRatios[ii], scaleRatios[ii]);
                    temp_validMaps = myResize(activeMap[t][k],scaleRatios[ii],scaleRatios[ii]);
                    dFResize[ii][t].emplace_back(temp_dF);
                    validMaps[ii][t].emplace_back(temp_validMaps);
                }//for(k)
            }//for(t)
            H0s[ii] = std::ceil(H/scaleRatios[ii]);
            W0s[ii] = std::ceil(W/scaleRatios[ii]);
        }//for(ii)

        //seed map
        std::vector<std::vector<cv::Mat>> zscoreMap(dF.size());
        for (int t = 0; t < dF.size(); ++t) {
            for (int k = 0; k < dF[0].size(); ++k) {
                zscoreMap[t].emplace_back(cv::Mat::zeros(dF[0][0].rows,dF[0][0].cols, CV_32F));
            }
        }//for(t)

        for (int ii = 0; ii < Thrs.size(); ++ii) { // threshold
            float curThr = Thrs[ii];
            for (int ii_ds = 0; ii_ds < scaleRatios.size(); ++ii_ds) { // downsample rate
                float H0 = H0s[ii_ds];
                float W0 = W0s[ii_ds];
                float scaleRatio = scaleRatios[ii_ds];
                std::vector<int> tmp;
                for (int ind = 0; ind < scaleRatio; ++ind) {
                    for (int ind_i = 0; ind_i < scaleRatio; ++ind_i) {
                        tmp.emplace_back(ind+1);
                    }
                }//for(ind)
                std::vector<std::vector<cv::Mat>> selectMap(dFResize[0].size());
                for (int t = 0; t < dFResize[0].size(); ++t) {
                    for (int k = 0; k < dFResize[0][0].size(); ++k) {
                        selectMap[t].emplace_back(cv::Mat::zeros(dFResize[0][0][0].rows, dFResize[0][0][0].cols, CV_32S));
                        for (int i = 0; i < dFResize[0][0][0].rows; ++i) {
                            for (int j = 0; j < dFResize[0][0][0].cols; ++j) {
                                if ((dFResize[ii_ds][t][k].at<float>(i,j) > curThr) && (validMaps[ii][t][k].at<int>(i,j) !=0)){
                                    selectMap[t][k].at<int>(i,j) = 1;
                                }//if
                            }//for(j)
                        }//for(i)
                    }//for(k)
                }//for
                std::vector<std::vector<Point_struct>> curRegions = bw2Reg(selectMap);
                // rough filter -- for acceleration
                for (int ii_cur = 0; ii_cur < curRegions.size(); ++ii_cur) {
                    if (curRegions[ii_cur].size() <= (opts.minSize / pow(scaleRatio,2) * opts.minDur / 3)){
                        curRegions.erase(curRegions.begin() + ii_cur);
                    }
                }
                for (int ii_cur = 0; ii_cur < curRegions.size(); ++ii_cur) {
                    std::unordered_set<int> ihw;
                    std::vector<int> ih;
                    std::vector<int> iw;
                    std::vector<int> il;
                    std::vector<int> it;
                    std::vector<Point_struct> pix;
                    int max_it = 0, min_it = 5000;
                    for (int ii_ite = 0; ii_ite < curRegions[ii_cur].size(); ++ii_ite) { //ii_ite -- different points in a region
                        Point_struct temp;
                        temp.i = curRegions[ii_cur][ii_ite].i;
                        temp.j = curRegions[ii_cur][ii_ite].j;
                        temp.k = curRegions[ii_cur][ii_ite].k;
                        temp.t = curRegions[ii_cur][ii_ite].t;
                        pix.emplace_back(temp);
                        ih.emplace_back(curRegions[ii_cur][ii_ite].i);
                        iw.emplace_back(curRegions[ii_cur][ii_ite].j);
                        il.emplace_back(curRegions[ii_cur][ii_ite].k);
                        it.emplace_back(curRegions[ii_cur][ii_ite].t);
                        int ihw_ind = sub2ind(ih[ii_ite],iw[ii_ite],il[ii_ite],H,W);
                        ihw.insert(ihw_ind);
                        }//for(ii_ite)
                    int dur = *std::max_element(ih.begin(),ih.end()) - *std::min_element(ih.begin(),ih.end()) + 1;
                    int arLabel_hs = ih[0] * scaleRatio;
                    int arLabel_he = std::min(static_cast<float>(H),ih[1]*scaleRatio - 1);
                    int arLabel_ws = iw[0] * scaleRatio;
                    int arLabel_we = std::min(static_cast<float>(W),iw[1]*scaleRatio - 1);
                    std::unordered_set<int> arLabel;
                    for (int i = arLabel_hs; i <= arLabel_he; ++i) {
                        for (int j = arLabel_ws; j < arLabel_we; ++j) {
                            arLabel.insert(activeMap[it[0]][il[0]].at<int>(i,j));
                        }
                    }
                    arLabel.erase(0);
                    int arLabel_val = *arLabel.begin();
                        //filter according to size and duration, also check seed detected or not
                    if (dur < opts.minDur || ihw.size()< std::max(static_cast<float>(opts.minSize), regSz[arLabel_val]*opts.seedSzRatio) / scaleRatio/scaleRatio){
                        continue;
                    }
                    // convert back
                    cv::Mat ihOrg = cv::Mat(ih.size(), scaleRatio*scaleRatio, CV_32F);
                    cv::Mat iwOrg = cv::Mat(ih.size(), scaleRatio*scaleRatio, CV_32F);
                    cv::Mat ilOrg = cv::Mat(ih.size(), scaleRatio*scaleRatio, CV_32F);
                    cv::Mat itOrg = cv::Mat(ih.size(), scaleRatio*scaleRatio, CV_32F);
                    cv::Mat select = cv::Mat::zeros(ih.size(), scaleRatio*scaleRatio, CV_8U);
                    for (int col = 0, add = 1; col < scaleRatio*scaleRatio; ++col,++add) {
                        if (add == scaleRatio+1){
                            add = 1;
                        }
                        for (int row = 0; row < ih.size(); ++row) {
                            ihOrg.at<float>(row,col) = ih[row]*scaleRatio + add;
                            iwOrg.at<float>(row,col) = iw[row]*scaleRatio + add;
                            ilOrg.at<float>(row,col) = il[row];
                            itOrg.at<float>(row,col) = it[row];
                        }
                    }//for(col)
                    for (int i = 0; i < ih.size(); ++i) {
                        for (int j = 0; j < scaleRatio * scaleRatio; ++j) {
                            if ((ihOrg.at<float>(i,j) <= H) && (iwOrg.at<float>(i,j) <= W)){
                                select.at<uint>(i,j) = 1;
                            }
                        }
                    }
                    std::vector<Point_struct> pixOrg;
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
                    if ((pixOrg.size()< std::max(static_cast<float>(opts.minSize), regSz[arLabel_val]*opts.seedSzRatio)) && find_flag){
                        continue;
                    }

                    //calculate significance
                    float t_scl =  std::max(static_cast<double>(1), std::round(dur/opts.TPatch));



                }//for(ii_cur)
            }//for(ii_ds) --scaleRatios.size()
        }//for(ii) --Thrs.size()
        




    }//seedDetect2_DS_accelerate


    void seDetection(const std::vector<std::vector<cv::Mat>>& dF, const std::vector<std::vector<cv::Mat>>& dataOrg,
                     const std::vector<std::vector<Point_struct>>& arLst){
        //may be according to active region size, set different scale
        for (int i = opts.maxSpaScale; i <= opts.minSpaScale; --i) {
            opts.scaleRatios.emplace_back(i);
        }
        std::cout<< "--------start seed detecting--------"<<std::endl;
        /*
         * test split by xuelong's methods from active region
         * top-down seed detection
         * is there any advantage that when checking seeds, calculate the score
         * for each pixel then combine them together? -- Since each individual
         * pixel in seed has different time windows. Consider the potential
         * propagation, we cannot use the average curve for seed detection.
         */




    }//seDetection()


    void phaseRun(){
        std::cout<< "--------start temporal segmentation--------"<<std::endl;
        if (opts.needTemp){
            opts.step = 0.5; // sigma
            opts.sigThr = 3.5; // Z score of seed significance
            opts.maxDelay = 0.6; // allowed maximum dissimilarity in merging
            opts.seedSzRatio = 0.01; // seed size relative to active region
            opts.needRefine = false; // peaks are temporally adjacent, need regine
            opts.needGrow = false; // grow active regions according to signal pattern

            std::vector<std::vector<cv::Mat>> dF1 = opts.dF1;
            std::vector<std::vector<cv::Mat>> dataOrg1 = opts.data1_org;
            std::vector<std::vector<Point_struct>> arLst1 = opts.arLst1;



        }//if

    }//phaseRun

}//namespace
