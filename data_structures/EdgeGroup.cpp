//
// Created by jure on 11. 08. 21.
//

#include "EdgeGroup.h"

void data_structures::EdgeGroup::addEdge(const std::shared_ptr<data_structures::Edge>& edge) {
    this->edges.push_back(edge);
}

std::vector<std::shared_ptr<data_structures::Edge>> data_structures::EdgeGroup::getEdges() {
    return this->edges;
}
