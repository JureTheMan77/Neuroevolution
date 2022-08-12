//
// Created by jure on 11. 03. 21.
//

#include "Agent.h"

std::shared_ptr<data_structures::Graph> evolution::Agent::getGraph() {
    return this->graph;
}

double evolution::Agent::getFitness() {
    return this->fitness;
}

void evolution::Agent::setFitness(double newFitness) {
    this->fitness = newFitness;
}

std::string evolution::Agent::toString(bool technical) {
    return "Fitness: " + std::to_string(this->fitness) + "\n" + this->graph->toString(technical);
}

std::shared_ptr<evolution::Agent>
evolution::Agent::create(const std::shared_ptr<data_structures::Graph> &graphArg) {
    return std::make_shared<evolution::Agent>(graphArg);
}

std::shared_ptr<evolution::Agent> evolution::Agent::deepClone() {
    // also copy fitness
    auto graph = this->getGraph()->deepClone();
    return std::make_shared<evolution::Agent>(graph, this->fitness, this->accuracy, this->matthewsCorrelationCoefficient);
}

std::shared_ptr<evolution::Agent>
evolution::Agent::create(unsigned int inputVertices, std::vector<std::string> &inputLabels, unsigned int outputVertices,
                         std::vector<std::string> &outputLabels, double maxMutationChance) {
    auto graph = data_structures::Graph::createGraph(inputVertices, inputLabels, 0, outputVertices, outputLabels);
    return std::make_shared<evolution::Agent>(graph);
}

bool evolution::Agent::isNewAgent() const {
    return newAgent;
}

void evolution::Agent::setNewAgent(bool newAgentArg) {
    Agent::newAgent = newAgentArg;
}

double evolution::Agent::getAccuracy() const {
    return accuracy;
}

void evolution::Agent::setAccuracy(double accuracy) {
    Agent::accuracy = accuracy;
}

double evolution::Agent::getMatthewsCorrelationCoefficient() const {
    return matthewsCorrelationCoefficient;
}

void evolution::Agent::setMatthewsCorrelationCoefficient(double matthewsCorrelationCoefficient) {
    Agent::matthewsCorrelationCoefficient = matthewsCorrelationCoefficient;
}
