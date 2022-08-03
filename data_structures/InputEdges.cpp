//
// Created by jure on 11. 08. 21.
//

#include <algorithm>
#include "InputEdges.h"

bool data_structures::InputEdges::allInputEdgesTraversed() {
    return std::all_of(this->edges.begin(), this->edges.end(),
                       [](const std::shared_ptr <data_structures::Edge>& edge) { return edge->isTraversed(); });
}

std::vector <std::shared_ptr<data_structures::Edge>> data_structures::InputEdges::getInputEdges() {
    return this->getEdges();
}

void data_structures::InputEdges::addInputEdge(const std::shared_ptr <data_structures::Edge> &edge) {
    this->addEdge(edge);
}
