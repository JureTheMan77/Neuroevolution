//
// Created by jure on 7/13/20.
//

#include <sstream>
#include <iomanip>
#include <limits>
#include "Edge.h"
#include "Vertex.h"

std::shared_ptr<data_structures::Vertex> data_structures::Edge::getInput() {
    return this->input.lock();
}

std::shared_ptr<data_structures::Vertex> data_structures::Edge::getOutput() {
    return this->output.lock();
}

std::shared_ptr<data_structures::Edge>
data_structures::Edge::createEdge(const std::shared_ptr<data_structures::Vertex> &inputArg,
                                  const std::shared_ptr<data_structures::Vertex> &outputArg,
                                  unsigned int indexArg,
                                  double weightArg,
                                  unsigned int traverseLimitArg, double mutationChance, bool dominant,
                                  double chanceToGetDominated, unsigned int maxChildren) {
    return std::make_shared<data_structures::Edge>(inputArg, outputArg, indexArg, weightArg, traverseLimitArg,
                                                   mutationChance, dominant, chanceToGetDominated, maxChildren);
}

double data_structures::Edge::propagateValue() {
    double outputValue = this->input.lock()->getValue() * this->weight;
    this->output.lock()->combineValue(outputValue);
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

    stream << std::setprecision(std::numeric_limits<double>::digits10) << "{" << this->input.lock()->toString() << " --{W:"
           << this->weight << " TL:" << this->traverseLimit << "}-> "
           << this->output.lock()->toString() << "}";

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

unsigned int data_structures::Edge::getIndex() const {
    return index;
}

void data_structures::Edge::setIndex(unsigned int index) {
    this->index = index;
}
