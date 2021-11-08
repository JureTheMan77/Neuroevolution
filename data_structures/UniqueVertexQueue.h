//
// Created by jure on 24. 02. 21.
//

#ifndef NEUROEVOLUTION_UNIQUEVERTEXQUEUE_H
#define NEUROEVOLUTION_UNIQUEVERTEXQUEUE_H

#include <queue>
#include <memory>
#include <unordered_set>
#include "Vertex.h"

namespace data_structures {
    class UniqueVertexQueue {
    private:
        std::queue<std::shared_ptr<data_structures::Vertex>> queue{};
        std::unordered_set<std::shared_ptr<data_structures::Vertex>> set{};
    public:
        /**
         * Creates a new empty queue object.
         */
        UniqueVertexQueue() = default;

        /**
         * Destroys this object.
         */
        ~UniqueVertexQueue() = default;

        /**
         * Adds @param element to the back of the queue if it has not been added before.
         * @param element element to add
         */
        void enqueue(const std::shared_ptr<data_structures::Vertex> &element);

        /**
         * Removes the first element from the queue.
         * @return the removed element
         */
        std::shared_ptr<data_structures::Vertex> dequeue();

        /**
         * Is the queue empty?
         * @return true, if the queue is empty
         */
        bool empty();

        /**
         * Removes all elements from this queue.
         */
        void clear();
    };
}


#endif //NEUROEVOLUTION_UNIQUEVERTEXQUEUE_H
