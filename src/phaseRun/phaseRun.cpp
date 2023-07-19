//
// Created by Kevin Qiao on 7/18/2023.
//

#include "phaseRun.h"

namespace AQuA{

    std::vector<std::vector<std::vector<std::vector<cv::Mat>>>> normalizeAndResize(const std::vector<std::vector<cv::Mat>>& dataOrg){
        std::vector<int> scaleRatios = opts.scaleRatios;
        std::vector<std::vector<cv::Mat>> datDS(T);
        std::vector<std::vector<std::vector<std::vector<cv::Mat>>>> datResize(scaleRatios.size());

        //rescl dForg
        for (int ii = 0; ii < scaleRatios.size(); ++ii) {
            int scaleRatio = scaleRatios[ii];
            for (int t = 0; t < T; ++t) {
                datDS[t].resize(L);
                for (int k = 0; k < L; ++k) {
                    cv::resize(dataOrg[t][k],datDS[t][k],cv::Size(),1/scaleRatio,1/scaleRatio, cv::INTER_AREA);
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
                cv::resize(opts.tempVarOrg1[k], var1[k], cv::Size(), 1 / scaleRatio, 1 / scaleRatio, cv::INTER_AREA);
                cv::resize(temp1, var2[k], cv::Size(), 1 / scaleRatio, 1 / scaleRatio, cv::INTER_AREA);

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
        std::vector<int> scaleRatios = {2,4,8};
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
            std::unordered_set<int> unique_indices;
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
        std::vector<std::vector<std::vector<std::vector<cv::Mat>>>> validMaps(scaleRatios.size()); //treat as boolean
        std::vector<std::vector<std::vector<std::vector<cv::Mat>>>> datResize = normalizeAndResize(dataOrg); //normalized data to do significance test
        std::vector<std::vector<std::vector<std::vector<cv::Mat>>>> dFResize(scaleRatios.size()); //down sampled data to do selection
        std::vector<int> H0s(scaleRatios.size(), 0);
        std::vector<int> W0s(scaleRatios.size(), 0);
        for (int ii = 0; ii < scaleRatios.size(); ++ii) {
//            datResize[ii] = reshape
            for (int t = 0; t < T; ++t) {
                for (int k = 0; k < L; ++k) {
                    cv::resize(dF[t][k], dFResize[ii][t][k], cv::Size(), 1 / scaleRatios[ii], 1 / scaleRatios[ii], cv::INTER_AREA);
                    cv::resize(activeMap[t][k], validMaps[ii][t][k], cv::Size(), 1 / scaleRatios[ii], 1 / scaleRatios[ii], cv::INTER_AREA);
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
