//
// Created by Kevin Qiao on 3/28/2023.
//


#include "flow.h"

namespace AQuA{

    int menu(){
        std::cout << "select the task: " << std::endl;
        std::cout << "1.preProcessRun" << std::endl;
        std::cout << "2.actRun" << std::endl;
        std::cout << "3.phaseRun" << std::endl;
        std::cout << "4.evtRun" << std::endl;
        std::cout << "5.feaRun" << std::endl;
        int ixTab = 1;
        std::cout << std::endl;
        std::cout << "input:  ";
        std::cin >> ixTab;
        return ixTab;
    }// menu()
}// namespace