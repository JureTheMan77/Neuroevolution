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

std::weak_ptr<data_structures::Vertex> data_structures::Edge::getWeakInput() {
    return this->input;
}

std::shared_ptr<data_structures::Vertex> data_structures::Edge::getOutput() {
    return this->output.lock();
}

std::weak_ptr<data_structures::Vertex> data_structures::Edge::getWeakOutput() {
    return this->output;
}

std::shared_ptr<data_structures::Edge>
data_structures::Edge::createEdge(const std::shared_ptr<data_structures::Vertex> &inputArg,
                                  const std::shared_ptr<data_structures::Vertex> &outputArg, unsigned int indexArg,
                                  double weightArg, unsigned int traverseLimitArg, double mutationChance, bool dominant,
                                  unsigned int maxChildren) {
    return std::make_shared<data_structures::Edge>(inputArg, outputArg, indexArg, weightArg, traverseLimitArg,
                                                   mutationChance, dominant, maxChildren);
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


std::string data_structures::Edge::toString(bool technical) {
    std::ostringstream stream;

    if (technical) {
        stream << std::setprecision(3) << this->input.lock()->toString(technical)
               << " "
               << this->output.lock()->toString(technical)
               << " \""
               << this->weight << "," << this->traverseLimit << "\"";
    } else {
        stream << std::setprecision(std::numeric_limits<double>::digits10) << "{"
               << this->input.lock()->toString(technical)
               << " --{W:"
               << this->weight << " TL:" << this->traverseLimit << "}-> "
               << this->output.lock()->toString(technical) << "}";
    }
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

bool data_structures::Edge::isFlaggedForDeletion() const {
    return flaggedForDeletion;
}

void data_structures::Edge::setFlaggedForDeletion(bool flaggedForDeletion) {
    Edge::flaggedForDeletion = flaggedForDeletion;
}

void data_structures::Edge::setWeight(double weightArg) {
    this->weight = weightArg;
}

void data_structures::Edge::setTraverseLimit(unsigned int traverseLimit) {
    this->traverseLimit = traverseLimit;
}
