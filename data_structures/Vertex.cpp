//
// Created by jure on 7/13/20.
//

#include "Vertex.h"

std::string data_structures::Vertex::toString() {
    return "";
}

//std::unique_ptr<data_structures::Vertex> data_structures::Vertex::createVertex(unsigned int vertexIndex) {
//    return std::make_unique<Vertex>(vertexIndex);
//}

void data_structures::Vertex::combineValue(double argValue) {
    this->value += argValue;
}

double data_structures::Vertex::getValue() {
    return this->value;
}

void data_structures::Vertex::reset() {
    this->value = 0;
}

unsigned int data_structures::Vertex::getIndex() {
    return this->index;
}

bool data_structures::Vertex::isDominant() {
    return this->dominant;
}

double data_structures::Vertex::getChanceToGetDominated() {
    return this->chanceToGetDominated;
}

unsigned int data_structures::Vertex::getMaxChildren() {
    return this->maxChildren;
}
