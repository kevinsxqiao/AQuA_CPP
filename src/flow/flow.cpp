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
        std::cout << "0.end" << std::endl;
        int ixTab = 1;
        while(1){
            std::cout << std::endl;
            std::cout << "input:  ";
            std::cin >> ixTab;
            switch (ixTab) {
                case 0: ixTab = 0;
                    break;
                case 1: ixTab = 1;
                    break;
                case 2: ixTab = 2;
                    break;
                case 3: ixTab = 3;
                    break;
                case 4: ixTab = 4;
                    break;
                case 5: ixTab = 5;
                    break;
                default:
                    std::cout << "Invalid input. Please enter a valid option." << std::endl;
                    continue;
            }
            return ixTab;
        }
    }// menu()

}// namespace