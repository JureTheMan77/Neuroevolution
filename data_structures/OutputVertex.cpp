//
// Created by jure on 4. 02. 21.
//

#include <sstream>
#include "OutputVertex.h"

std::string data_structures::OutputVertex::toString() {
    std::ostringstream stream;
    stream << "{Index: " << this->getIndex()
           << ", Label: " << this->label
           << ", Value: " << this->getValue()
           << "}";

    return stream.str();
}

std::shared_ptr<data_structures::OutputVertex>
data_structures::OutputVertex::createOutputVertex(unsigned int index, std::string label) {
    return std::make_shared<data_structures::OutputVertex>(index, label);
}

void data_structures::OutputVertex::addOutputEdge(const std::shared_ptr<data_structures::Edge> &edge) {
    // ignored
}
