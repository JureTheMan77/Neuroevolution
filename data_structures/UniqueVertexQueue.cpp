//
// Created by jure on 24. 02. 21.
//

#include "UniqueVertexQueue.h"

void data_structures::UniqueVertexQueue::enqueue(const std::shared_ptr <data_structures::Vertex> &element) {
    if (this->set.find(element) == this->set.end()) {
        this->queue.push(element);
        this->set.insert(element);
    }
}

std::shared_ptr <data_structures::Vertex> data_structures::UniqueVertexQueue::dequeue() {
    std::shared_ptr <data_structures::Vertex> element = this->queue.front();
    this->queue.pop();
    this->set.erase(element);
    return element;
}

bool data_structures::UniqueVertexQueue::empty() {
    return this->queue.empty();
}

void data_structures::UniqueVertexQueue::clear() {
    std::queue < std::shared_ptr < data_structures::Vertex >> ().swap(this->queue);
    this->set.clear();
}
