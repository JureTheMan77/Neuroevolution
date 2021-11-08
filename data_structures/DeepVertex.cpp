//
// Created by jure on 4. 02. 21.
//

#include <sstream>
#include "DeepVertex.h"

std::shared_ptr<data_structures::DeepVertex>
data_structures::DeepVertex::createDeepVertex(unsigned int index, bool dominant, double chanceToGetDominated,
                                              unsigned int maxChildren) {
    return std::make_shared<data_structures::DeepVertex>(index, dominant, chanceToGetDominated, maxChildren);
}

std::string data_structures::DeepVertex::toString() {
    std::ostringstream stream;
    stream << "{Index: " << this->getIndex()
           << ", Value: " << this->getValue()
           << "}";

    return stream.str();
}


