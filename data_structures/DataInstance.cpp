//
// Created by jure on 11. 02. 21.
//

#include <vector>
#include <memory>
#include "DataInstance.h"

std::vector<double> data_structures::DataInstance::getValues() {
    return this->values;
}

unsigned int data_structures::DataInstance::getCorrectIndex() {
    return this->correctIndex;
}

std::shared_ptr<data_structures::DataInstance>
data_structures::DataInstance::createDataInstance(const std::vector<double>& valuesAndLabelIndex) {
    return std::make_shared<data_structures::DataInstance>(valuesAndLabelIndex);
}

