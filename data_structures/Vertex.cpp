//
// Created by jure on 7/13/20.
//

#include "Vertex.h"

std::string data_structures::Vertex::toString([[maybe_unused]] bool technical) {
    return "";
}

void data_structures::Vertex::combineValue(double argValue) {
    this->value += argValue;
    this->updateActivationValue();
}

double data_structures::Vertex::getValue() const {
    return this->value;
}

void data_structures::Vertex::reset() {
    this->value = 0;
    this->setVisited(false);
}

unsigned int data_structures::Vertex::getIndex() {
    return this->index;
}

bool data_structures::Vertex::isVisited() const {
    return visited;
}

void data_structures::Vertex::setVisited(bool visitedArg) {
    this->visited = visitedArg;
}

void data_structures::Vertex::setIndex(unsigned int indexArg) {
    this->index = indexArg;
}

void data_structures::Vertex::updateActivationValue() {
    this->activationValue = this->value;
}

double data_structures::Vertex::getActivationValue() const {
    return this->activationValue;
}
