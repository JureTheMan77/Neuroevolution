//
// Created by jure on 3. 02. 21.
//

#include <sstream>
#include "InputVertex.h"

std::string data_structures::InputVertex::toString() {
    std::ostringstream stream;
    stream << "{Index: " << this->getIndex()
           << ", Label: " << this->label
           << ", Value: " << this->getValue()
           << "}";

    return stream.str();
}

std::shared_ptr<data_structures::InputVertex>
data_structures::InputVertex::createInputVertex(unsigned int index, bool dominant, double chanceToGetDominated,
                                                unsigned int maxChildren, std::string &labelArg) {
    return std::make_shared<data_structures::InputVertex>(index, dominant, chanceToGetDominated, maxChildren, labelArg);
}

void data_structures::InputVertex::addInputEdge(const std::shared_ptr<data_structures::Edge> &edge) {
    // ignored
}

bool data_structures::InputVertex::allInputEdgesTraversed() {
    return true;
}
