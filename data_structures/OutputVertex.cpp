//
// Created by jure on 4. 02. 21.
//

#include <sstream>
#include "OutputVertex.h"

std::string data_structures::OutputVertex::toString(bool technical) {
    std::ostringstream stream;
    if (technical) {
        stream << "\"" << this->label << "\"";
    } else {
        stream << "{Index: " << this->getIndex()
               << ", Label: " << this->label
               //<< ", Value: " << this->getValue()
               << "}";
    }

    return stream.str();
}

std::shared_ptr<data_structures::OutputVertex>
data_structures::OutputVertex::createOutputVertex(unsigned int index, std::string label) {
    return std::make_shared<data_structures::OutputVertex>(index, label);
}

void data_structures::OutputVertex::addOutputEdge(const std::shared_ptr<data_structures::Edge> &edge) {
    // ignored
}

std::shared_ptr<data_structures::OutputVertex> data_structures::OutputVertex::deepClone() {
    return std::make_shared<data_structures::OutputVertex>(this->getIndex(), this->label);
}

enums::VertexType data_structures::OutputVertex::getType() {
    return enums::VertexType::Output;
}

const std::string &data_structures::OutputVertex::getLabel() const {
    return label;
}
