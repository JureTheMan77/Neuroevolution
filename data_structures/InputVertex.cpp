//
// Created by jure on 3. 02. 21.
//

#include <sstream>
#include "InputVertex.h"

std::string data_structures::InputVertex::toString(bool technical) {
    std::ostringstream stream;
    if (technical) {
        stream << "\"" << this->label << "\"";
    } else {
        stream << "{Index: " << this->getIndex()
               << ", Label: " << this->label
               << ", Value: " << this->getActivationValue()
               << "}";
    }

    return stream.str();
}

std::shared_ptr<data_structures::InputVertex>
data_structures::InputVertex::createInputVertex(unsigned int index, std::string &labelArg) {
    return std::make_shared<data_structures::InputVertex>(index, labelArg);
}

void data_structures::InputVertex::addInputEdge(const std::shared_ptr<data_structures::Edge> &edge) {
    // ignored
}

bool data_structures::InputVertex::allInputEdgesTraversedOnce() {
    return true;
}

std::shared_ptr<data_structures::InputVertex> data_structures::InputVertex::deepClone() {
    return std::make_shared<data_structures::InputVertex>(this->getIndex(), this->label);
}

enums::VertexType data_structures::InputVertex::getType() {
    return enums::VertexType::Input;
}

const std::string &data_structures::InputVertex::getLabel() const {
    return label;
}
