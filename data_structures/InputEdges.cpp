//
// Created by jure on 11. 08. 21.
//

#include <algorithm>
#include <iostream>
#include "InputEdges.h"

bool data_structures::InputEdges::allInputEdgesTraversedOnce() {
    //return std::all_of(this->edges.begin(), this->edges.end(),
    //                   [](const std::shared_ptr<data_structures::Edge> &edge) { return edge->isTraversedOnce(); });
    bool result = true;
    for (const auto &edge: this->edges) {
        if (!edge->isTraversedOnce()) {
            result = false;
            break;
        }
    }
    return result;
}

std::vector<std::shared_ptr<data_structures::Edge>> data_structures::InputEdges::getInputEdges() {
    return this->getEdges();
}

void data_structures::InputEdges::addInputEdge(const std::shared_ptr<data_structures::Edge> &edge) {
    if(edge == nullptr){
        std::cout<<"InputEdges.cpp:27";
        std::invalid_argument("Edge is null!");
    }
    this->addEdge(edge);
}

void data_structures::InputEdges::eraseInputEdge(unsigned long position) {
    this->eraseEdge(position);
}

unsigned int data_structures::InputEdges::inputEdgeTraversalsRemaining() {
    unsigned int sum = 0;
    for (const auto &edge: this->edges) {
        sum += edge->traversalsRemaining();
    }
    return sum;
}
