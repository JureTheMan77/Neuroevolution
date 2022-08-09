//
// Created by jure on 11. 08. 21.
//

#include "EdgeGroup.h"

void data_structures::EdgeGroup::addEdge(const std::shared_ptr<data_structures::Edge> &edge) {
    this->edges.push_back(edge);
}

std::vector<std::shared_ptr<data_structures::Edge>> data_structures::EdgeGroup::getEdges() {
    return this->edges;
}

void data_structures::EdgeGroup::replaceEdges(std::vector<std::shared_ptr<data_structures::Edge>> const &newEdges) {
    this->edges.clear();
    for(const auto &edge : newEdges){
        this->edges.push_back(edge);
    }
}

void data_structures::EdgeGroup::clearEdges() {
    this->edges.clear();
}

void data_structures::EdgeGroup::eraseEdge(unsigned long position) {
    this->edges.erase(this->edges.begin() + position);
}
