//
// Created by jure on 4. 02. 21.
//

#include <sstream>
#include "DeepVertex.h"

std::shared_ptr<data_structures::DeepVertex>
data_structures::DeepVertex::createDeepVertex(unsigned int index, bool dominant, double chanceToGetDominated,
                                              double mutationChanceArg, unsigned int maxChildren) {
    return std::make_shared<data_structures::DeepVertex>(index, dominant, chanceToGetDominated, mutationChanceArg,
                                                         maxChildren);
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
    return std::make_shared<data_structures::DeepVertex>(this->getIndex(), this->isDominant(),
                                                         this->getChanceToGetDominated(),
                                                         this->getMutationChance(),
                                                         this->getMaxChildren());
}

enums::VertexType data_structures::DeepVertex::getType() {
    return enums::VertexType::Deep;
}

bool data_structures::DeepVertex::isFlaggedForDeletion() const {
    return flaggedForDeletion;
}

void data_structures::DeepVertex::setFlaggedForDeletion(bool flaggedForDeletion) {
    DeepVertex::flaggedForDeletion = flaggedForDeletion;
}

bool data_structures::DeepVertex::allInputEdgesFlaggedForDeletion() {
    for (const auto &edge: this->getInputEdges()) {
        if (!edge->isFlaggedForDeletion()) {
            return false;
        }
    }
    return true;
}

bool data_structures::DeepVertex::allOutputEdgesFlaggedForDeletion() {
    for (const auto &edge: this->getOutputEdges()) {
        if (!edge->isFlaggedForDeletion()) {
            return false;
        }
    }
    return true;
}

void data_structures::DeepVertex::combineValue(double argValue) {
    // leaky RELU
    double computedValue = argValue;
    if (argValue < 0) {
        computedValue = computedValue * 0.01;
    }

    Vertex::combineValue(computedValue);
}
