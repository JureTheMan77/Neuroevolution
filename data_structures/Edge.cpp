//
// Created by jure on 7/13/20.
//

#include <sstream>
#include "Edge.h"
#include "Vertex.h"

std::shared_ptr<data_structures::Vertex> data_structures::Edge::getInput() {
    return this->input;
}

std::shared_ptr<data_structures::Vertex> data_structures::Edge::getOutput() {
    return this->output;
}

std::shared_ptr<data_structures::Edge>
data_structures::Edge::createEdge(std::shared_ptr<data_structures::Vertex> inputArg,
                                  std::shared_ptr<data_structures::Vertex> outputArg, double weightArg,
                                  unsigned int traverseLimitArg) {
    return std::make_shared<data_structures::Edge>(inputArg, outputArg, weightArg, traverseLimitArg);
}

double data_structures::Edge::propagateValue() {
    double outputValue = this->input->getValue() * this->weight;
    this->output->combineValue(outputValue);
    this->traverseCount += 1;
    return outputValue;
}

bool data_structures::Edge::isTraversed() {
    return this->traverseCount > 0;
}

bool data_structures::Edge::isAtTaversalLimit() {
    return this->traverseLimit == this->traverseCount;
}


std::string data_structures::Edge::toString() {
    std::ostringstream stream;

    stream << "{" << this->input->toString() << " --{W:" << this->weight << " TL:" << this->traverseLimit << "}-> "
           << this->output->toString() << "}";

    return stream.str();
}

void data_structures::Edge::reset() {
    this->traverseCount = 0;
}

void data_structures::Edge::combineWeight(double argWeight) {
    this->weight += argWeight;
}

double data_structures::Edge::getWeight() {
    return weight;
}

unsigned int data_structures::Edge::getTraverseLimit() {
    return traverseLimit;
}
