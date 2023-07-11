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
        std::cout<<"connected regions: "<<regions.size()<<std::endl;
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


    void acDetect(std::vector<std::vector<cv::Mat>> dF1, bool*** evtSpatialMask) {
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
        std::vector<std::vector<cv::Mat>> activeMap(T);
        std::vector<std::vector<cv::Mat>> selectMap(T);
        for (int t = 0; t < T; ++t) {
            for (int k = 0; k < L; ++k) {
                activeMap[t][k] = cv::Mat::zeros(H,W,CV_32S);
                selectMap[t][k] = cv::Mat::zeros(H,W,CV_32S);
            }
        }
        float nReg = 0;
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
        }//for(k_thr)

    }



    /*
     * active region detection and update overlay map
     */
    void actRun(){
        std::cout<< "--------start detecting--------"<<std::endl;
        bool*** evtSpatialMask = createEvtSpatialMask();

        //foreground and seed detection
        acDetect(opts.dF1, evtSpatialMask);

    }

}
