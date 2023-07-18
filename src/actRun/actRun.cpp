//
// Created by Kevin Qiao on 7/6/2023.
//
#include "actRun.h"

namespace AQuA{

    /*
     * std::vector<std::vector<std::vector<int>>> regions;
     * the innermost Point_struct is 4d point location = {t,k,i,j};
     * the outer vector contains different points
     * the outermost vector contains different regions
     */
    std::vector<std::vector<Point_struct>> bwconncomp4D(std::vector<std::vector<cv::Mat>> BW){
        std::vector<std::vector<cv::Mat>> visited(BW.size(),std::vector<cv::Mat>(BW[0].size()));
        for (int t = 0; t < BW.size(); ++t) {
            for (int k = 0; k < BW[0].size(); ++k) {
                visited[t][k] = cv::Mat::zeros(BW[0][0].rows, BW[0][0].cols, CV_32S);
            }
        }
        std::vector<std::vector<Point_struct>> regions;
        int dx[3] = {-1, 0, 1};
        int dy[3] = {-1, 0, 1};
        int dz[3] = {-1, 0, 1};
        int dt[3] = {-1, 0, 1};
        for (int t = 0; t < BW.size(); ++t) {
            for (int k = 0; k < BW[0].size(); ++k) {
                for (int i = 0; i < BW[0][0].rows; ++i) {
                    for (int j = 0; j < BW[0][0].cols; ++j) {
                        if ((BW[t][k].at<int>(i,j) == 1) && (visited[t][k].at<int>(i,j) != 1)){
                            std::vector<Point_struct> component;
                            Point_struct point;
                            point = {t,k,i,j};
                            std::queue<Point_struct> queue;
                            queue.push(point);
                            visited[t][k].at<int>(i,j) = 1;
                            while (!queue.empty()){
                                Point_struct current;
                                current = queue.front();
                                queue.pop();
                                component.push_back(current);
                                for (int dti = 0; dti < 3; ++dti) {
                                    for (int dki = 0; dki < 3; ++dki) {
                                        for (int dxi = 0; dxi < 3; ++dxi) {
                                            for (int dyi = 0; dyi < 3; ++dyi) {
                                                int nt = current.t + dt[dti];
                                                int nz = current.k + dz[dki];
                                                int nx = current.i + dx[dxi];
                                                int ny = current.j + dy[dyi];

                                                if ((nx >= 0) && (nx < BW[0][0].rows) && (ny >= 0) && (ny < BW[0][0].cols) && (nz >= 0) &&
                                                (nz < BW[0].size()) && (nt >= 0) && (nt < BW.size()) &&
                                                (BW[nt][nz].at<int>(nx,ny)==1) && (visited[nt][nz].at<int>(nx,ny)!=1)) {

                                                        Point_struct curr;
                                                        curr = {nt,nz,nx,ny};
                                                        queue.push(curr);
                                                        visited[nt][nz].at<int>(nx,ny) = 1;
                                                }//if
                                            }//for(dyi)
                                        }
                                    }
                                }
                            }//while
                            regions.push_back(component);
                        }//if
                    }//for(j)
                }
            }
        }//for(t)
//        std::cout<<"connected regions: "<<regions.size()<<std::endl;
        return regions;
    }//bwconncomp4D


    std::vector<std::vector<Point_struct>> bw2Reg(std::vector<std::vector<cv::Mat>> BW){
//        if (opts.spaMergeDist>0){
//            if (BW[0].size() == 1){
//
//            } else{
//
//            }
//            region
//            for (int i = 0; i < ; ++i) {
//
//            }
//            sz
//            regions
//        } else{
//
//        }
        return bwconncomp4D(BW);

    }


    std::vector<std::vector<Point_struct>> acDetect(std::vector<std::vector<cv::Mat>> dF1, bool*** evtSpatialMask) {
//        if (ch == 1){
//
//        }
        std::vector<float> thrs;
        if ((opts.thrARScl > opts.maxdF1) || (opts.maxSize >= H * W * L) && (opts.circularityThr == 0)) {
            //no advanced filter setting, single threshold
            thrs.push_back(opts.thrARScl);
        } else {
            //have advanced filter setting, multiple threshold
            float step = (opts.maxdF1 - opts.thrARScl) / 10;
            for (float thr = opts.thrARScl; thr <= opts.maxdF1; thr += step) {
                thrs.push_back(thr);
            }
        }//else
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                for (int k = 0; k < L; ++k) {
                    if (!evtSpatialMask[i][j][k]){
                        for (int t = 0; t < T; ++t) {
                            dF1[t][k].at<float>(i,j) = -1;
                        }//for(t)
                    }//if
                }//for(k)
            }
        }//for(i)

        //valid region
        std::vector<std::vector<cv::Mat>> activeMap(T, std::vector<cv::Mat>(L));
        std::vector<std::vector<cv::Mat>> selectMap(T, std::vector<cv::Mat>(L));
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                activeMap[t][k] = cv::Mat::zeros(H,W,CV_32S);
                selectMap[t][k] = cv::Mat::zeros(H,W,CV_32S);
            }
        }
        int nReg = 0;
        std::vector<std::vector<Point_struct>> arLst;
        for (int k_thr = 0; k_thr < thrs.size(); ++k_thr) {
            float thr = thrs[k_thr];
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    for (int i = 0; i < H; ++i) {
                        for (int j = 0; j < W; ++j) {
                            if ((dF1[t][k].at<float>(i,j) > thr) && (activeMap[t][k].at<int>(i,j) == 0)){
                                selectMap[t][k].at<int>(i,j) = 1;
                            } else{
                                selectMap[t][k].at<int>(i,j) = 0;
                            }
                        }
                    }
                }
            }//for(t)
            std::vector<std::vector<Point_struct>> curRegions;
            curRegions = bw2Reg(selectMap);
            std::vector<bool> valid(curRegions.size(), false);
            for (int i_cur = 0; i_cur < curRegions.size(); ++i_cur) {
                std::vector<int> ih;
                std::vector<int> iw;
                std::vector<int> il;
                std::vector<int> it;
                for (int i_ite = 0; i_ite < curRegions[i_cur].size(); ++i_ite) {
                    ih.push_back(curRegions[i_cur][i_ite].i);
                    iw.push_back(curRegions[i_cur][i_ite].j);
                    il.push_back(curRegions[i_cur][i_ite].k);
                    it.push_back(curRegions[i_cur][i_ite].t);
                }//for(i_ite)
//                int curSz = 0;
//                for (int i_ite = 0; i_ite < curRegions[i_cur].size(); ++i_ite) {
//                    if (curRegions[i_cur].size() == 1){
//                        if (curRegions[i_cur][i_ite])
//                    }
//                }
                int ih_min = *std::min_element(ih.begin(),ih.end());
                int H0 = *std::max_element(ih.begin(),ih.end()) - ih_min + 1;
                int iw_min = *std::min_element(iw.begin(),iw.end());
                int W0 = *std::max_element(iw.begin(),iw.end()) - iw_min + 1;
                int il_min = *std::min_element(il.begin(),il.end());
                int L0 = *std::max_element(il.begin(),il.end()) - il_min + 1;
                int it_min = *std::min_element(it.begin(),it.end());
                int T0 = *std::max_element(it.begin(),it.end()) - it_min + 1;
                for (int i_ite = 0; i_ite < curRegions[i_cur].size(); ++i_ite) {
                    ih[i_ite] = ih[i_ite] - ih_min;
                    iw[i_ite] = iw[i_ite] - iw_min;
                    il[i_ite] = il[i_ite] - il_min;
                    it[i_ite] = it[i_ite] - it_min;
                }//for(i_ite)
                int**** curMap = create4dMatrix_int(H0,W0,L0,T0);  //initialize all with 0
                int*** curMap_new = create3dMatrix_int(H0,W0,L0);//initialize all with 0
                int curSz = 0;
                for (int i_ite = 0; i_ite < curRegions[i_cur].size(); ++i_ite) {
                    curMap[ih[i_ite]][iw[i_ite]][il[i_ite]][it[i_ite]] = 1;
                }//for(i_ite)
                for (int i = 0; i < H0; ++i) {
                    for (int j = 0; j < W0; ++j) {
                        for (int k = 0; k < L0; ++k) {
                            int sum = 0;
                            for (int t = 0; t < T0; ++t) {
                                sum+= curMap[i][j][k][t];
                            }//for(t)
                            if (sum>=opts.compress*T0){
                                curMap_new[i][j][k] = 1;
                                ++curSz;
                            }//if
                        }//for(k)
                    }
                }//for(i)
                if (curSz > opts.maxSize || curSz < opts.minSize || T0 < opts.minDur){
                    continue;
                }
                if (opts.circularityThr == 0){
                    valid[i_cur] = true;
                    continue;
                }
                if (dF1[0].size() == 1){
//                    erodeMap = imerode(curMap,strel('disk',1));
//                    boundary = curSz - sum(erodeMap(:));
//                    circularity = 4*pi*curSz/(boundary^2);
                } else{

                }

//                if (circularity > opts.circularityThr){
//
//                }

            }//for(i_cur)

            int valid_count = 0;
            for (int i_cur = 0; i_cur < curRegions.size(); ++i_cur) {// different regions
                if (valid[i_cur]){
                    ++valid_count;
                    std::vector<Point_struct> region_result;
                    for (int i_ite = 0; i_ite < curRegions[i_cur].size(); ++i_ite) {// different points in a region
                        int tt = curRegions[i_cur][i_ite].t;
                        int kk = curRegions[i_cur][i_ite].k;
                        int ii = curRegions[i_cur][i_ite].i;
                        int jj = curRegions[i_cur][i_ite].j;
                        activeMap[tt][kk].at<int>(ii,jj) = nReg + i_cur + 1;
                        Point_struct point;
                        point.t = tt;
                        point.k = kk;
                        point.i = ii;
                        point.j = jj;
                        region_result.push_back(point);
                    }//for(i_ite)
                    arLst.push_back(region_result);
                }//if
            }//for(i_cur)
            nReg += valid_count;
        }//for(k_thr)
//        std::cout<<"number of regions detected: "<<arLst.size()<<std::endl;
        opts.arLst1 = arLst;
        return arLst;
    }//acDetect()



    /*
     * active region detection and update overlay map
     */
    void actRun(){
        std::cout<< "--------start active region detecting--------"<<std::endl;
        bool*** evtSpatialMask = createEvtSpatialMask();
//        std::vector<std::vector<cv::Mat>> dF1 = load4DData_clean("C:/Users/Kevin Qiao/Desktop/AQuA_data/dF_real.mat","dF");
        //foreground and seed detection
        acDetect(opts.dF1, evtSpatialMask);
        release3dMatrix_bool(evtSpatialMask,H,W);
    }

}
