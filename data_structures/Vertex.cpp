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

void data_structures::Vertex::setValue(double value) {
    this->value = value;
}

bool data_structures::Vertex::isVisited() const {
    return visited;
}

void data_structures::Vertex::setVisited(bool visited) {
    this->visited = visited;
}

void data_structures::Vertex::setIndex(unsigned int index) {
    this->index = index;
}
