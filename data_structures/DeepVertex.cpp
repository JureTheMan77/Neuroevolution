//
// Created by jure on 4. 02. 21.
//

#include <sstream>
#include <algorithm>
#include "DeepVertex.h"

std::shared_ptr<data_structures::DeepVertex>
data_structures::DeepVertex::createDeepVertex(unsigned int index, bool dominant) {
    return std::make_shared<data_structures::DeepVertex>(index, dominant);
}

std::string data_structures::DeepVertex::toString(bool technical) {
    std::ostringstream stream;
    if (technical) {
        stream << this->getIndex();
    } else {
        stream << "{Index: " << this->getIndex()
               //<< ", Value: " << this->getValue()
               << "}";
    }

    return stream.str();
}

std::shared_ptr<data_structures::DeepVertex> data_structures::DeepVertex::deepClone() {
    return std::make_shared<data_structures::DeepVertex>(this->getIndex(), this->isDominant());
}

enums::VertexType data_structures::DeepVertex::getType() {
    return enums::VertexType::Deep;
}

bool data_structures::DeepVertex::isFlaggedForDeletion() const {
    return this->flaggedForDeletion;
}

void data_structures::DeepVertex::setFlaggedForDeletion(bool flaggedForDeletionArg) {
    this->flaggedForDeletion = flaggedForDeletionArg;
}

bool data_structures::DeepVertex::allInputEdgesFlaggedForDeletion() {
    return std::all_of(this->getInputEdges().begin(), this->getInputEdges().end(),
                       [](const std::shared_ptr<data_structures::Edge> &edge) { return edge->isFlaggedForDeletion(); });
}

bool data_structures::DeepVertex::allOutputEdgesFlaggedForDeletion() {
    return std::all_of(this->getOutputEdges().begin(), this->getOutputEdges().end(),
                       [](const std::shared_ptr<data_structures::Edge> &edge) { return edge->isFlaggedForDeletion(); });
}

void data_structures::DeepVertex::combineValue(double argValue) {
    // leaky RELU
    double computedValue;
    if (argValue < 0) {
        computedValue = argValue * 0.01;
    } else {
        computedValue = argValue;
    }

    Vertex::combineValue(computedValue);
}
