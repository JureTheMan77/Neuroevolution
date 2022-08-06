//
// Created by jure on 24. 01. 22.
//

#include "ICrossoverable.h"

double data_structures::ICrossoverable::getMutationChance() const {
    return mutationChance;
}

bool data_structures::ICrossoverable::isDominant() const {
    return dominant;
}

unsigned int data_structures::ICrossoverable::getMaxChildren() const {
    return maxChildren;
}

void data_structures::ICrossoverable::setDominant(bool dominant) {
    this->dominant = dominant;
}
