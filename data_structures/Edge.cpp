//
// Created by jure on 7/13/20.
//

#include <sstream>
#include <iomanip>
#include <limits>
#include "Edge.h"
#include "Vertex.h"

std::shared_ptr<data_structures::Vertex> data_structures::Edge::getInput() const {
    return this->input.lock();
}

std::shared_ptr<data_structures::Vertex> data_structures::Edge::getOutput() const {
    return this->output.lock();
}

std::shared_ptr<data_structures::Edge>
data_structures::Edge::createEdge(const std::shared_ptr<data_structures::Vertex> &inputArg,
                                  const std::shared_ptr<data_structures::Vertex> &outputArg, unsigned int indexArg,
                                  double weightArg, unsigned int traverseLimitArg, bool dominant) {
    return std::make_shared<data_structures::Edge>(inputArg, outputArg, indexArg, weightArg, traverseLimitArg,
                                                   dominant);
}

double data_structures::Edge::propagateValue() {
    if (this->isAtTraverseLimit()) {
        return 0;
    }
    double outputValue = this->input.lock()->getActivationValue() * this->weight;
    this->output.lock()->combineValue(outputValue);
    this->traverseCount += 1;
    return outputValue;
}

bool data_structures::Edge::isTraversedOnce() const {
    return this->traverseCount > 0;
}

bool data_structures::Edge::isAtTraverseLimit() const {
    return this->traverseLimit <= this->traverseCount;
}


std::string data_structures::Edge::toString(bool technical) const {
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

double data_structures::Edge::getWeight() const {
    return weight;
}

unsigned int data_structures::Edge::getTraverseLimit() const {
    return traverseLimit;
}

unsigned int data_structures::Edge::getIndex() const {
    return index;
}

void data_structures::Edge::setIndex(unsigned int indexArg) {
    this->index = indexArg;
}

bool data_structures::Edge::isFlaggedForDeletion() const {
    return flaggedForDeletion;
}

void data_structures::Edge::setFlaggedForDeletion(bool flaggedForDeletionArg) {
    Edge::flaggedForDeletion = flaggedForDeletionArg;
}

void data_structures::Edge::setWeight(double weightArg) {
    this->weight = weightArg;
}

void data_structures::Edge::setTraverseLimit(unsigned int traverseLimitArg) {
    this->traverseLimit = traverseLimitArg;
}

unsigned int data_structures::Edge::traversalsRemaining() const {
    if (this->isAtTraverseLimit()) {
        return 0;
    } else {
        return this->traverseLimit - this->traverseCount;
    }
}
