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

double data_structures::ICrossoverable::getChanceToGetDominated() const {
    return chanceToGetDominated;
}

unsigned int data_structures::ICrossoverable::getMaxChildren() const {
    return maxChildren;
}