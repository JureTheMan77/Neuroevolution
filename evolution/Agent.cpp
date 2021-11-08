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

std::string evolution::Agent::toString() {
    return "Fitness: " + std::to_string(this->fitness) + "\n" + this->graph->toString();
}

std::shared_ptr<evolution::Agent>
evolution::Agent::create(const std::shared_ptr<data_structures::Graph> &graphArg, double fitnessArg) {
    return std::make_shared<evolution::Agent>(graphArg, fitnessArg);
}

std::shared_ptr<evolution::Agent>
evolution::Agent::create(const std::shared_ptr<data_structures::Graph> &graphArg) {
    return std::make_shared<evolution::Agent>(graphArg);
}

std::shared_ptr<evolution::Agent> evolution::Agent::copy(std::shared_ptr<Agent> agent) {
    return std::make_shared<evolution::Agent>(agent);
}
