//
// Created by jure on 18. 08. 21.
//

#include "OutputEdges.h"

void data_structures::OutputEdges::addOutputEdge(const std::shared_ptr<data_structures::Edge> &edge) {
    this->addEdge(edge);
}

std::vector<std::shared_ptr<data_structures::Edge>> data_structures::OutputEdges::getOutputEdges() {
    return this->getEdges();
}

void data_structures::OutputEdges::eraseOutputEdge(unsigned long position) {
    this->eraseEdge(position);
}

bool data_structures::OutputEdges::canTraverseOut() {
    for (const auto &edge: this->edges) {
        if (!edge->isAtTraverseLimit()) {
            return true;
        }
    }
    return false;
}
