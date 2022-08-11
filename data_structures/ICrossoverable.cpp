//
// Created by jure on 24. 01. 22.
//

#include "ICrossoverable.h"

bool data_structures::ICrossoverable::isDominant() const {
    return dominant;
}

void data_structures::ICrossoverable::setDominant(bool dominant) {
    this->dominant = dominant;
}
