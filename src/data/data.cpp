//
// Created by Kevin Qiao on 3/28/2023.
//

#include "data.h"

namespace AQuA{

    void getData(){

    }

    bool isDefault() { //judge if registration and bleach have been executed; ---- true = no ; false = both executed
        return (preSetting::registrateCorrect == preSetting::registrateCorrect_default
                && preSetting::bleachCorrect == preSetting::bleachCorrect_default);
    }

    void preProcessRun(){

    }

}// namespace