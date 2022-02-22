//
// Created by jure on 7/13/20.
//

#ifndef NEUROEVOLUTION_EDGE_H
#define NEUROEVOLUTION_EDGE_H

#include <memory>
#include "Edge.fwd.h"
#include "Vertex.fwd.h"
#include "ICrossoverable.h"

namespace data_structures {
    class Edge : public data_structures::ICrossoverable {
    private:
        std::weak_ptr<data_structures::Vertex> input;
        std::weak_ptr<data_structures::Vertex> output;
        unsigned int index{};
        double weight;
        unsigned int traverseLimit;
        unsigned int traverseCount{0};
        bool flaggedForDeletion = false;
    public:
        /**
         * Creates a new edge.
         * @param input input vertex
         * @param output output vertex
         * @param weight weight of the edge
         * @param traverseLimit how many times can this edge be traversed, minimum is 1
         */
        Edge(std::shared_ptr<data_structures::Vertex> input, std::shared_ptr<data_structures::Vertex> output,
             unsigned int index, double weight, unsigned int traverseLimit, double mutationChance, bool dominant,
             double chanceToGetDominated, unsigned int maxChildren) : ICrossoverable(mutationChance, dominant,
                                                                                     chanceToGetDominated, maxChildren),
                                                                      input(std::move(input)),
                                                                      output(std::move(output)),
                                                                      index(index),
                                                                      weight(weight),
                                                                      traverseLimit(traverseLimit) {};

        static std::shared_ptr<data_structures::Edge> createEdge(const std::shared_ptr<data_structures::Vertex> &input,
                                                                 const std::shared_ptr<data_structures::Vertex> &output,
                                                                 unsigned int index,
                                                                 double weight, unsigned int traverseLimit,
                                                                 double mutationChance, bool dominant,
                                                                 double chanceToGetDominated, unsigned int maxChildren);

        /**
         * Destroys this object.
         */
        ~Edge() = default;

        /**
         * Get the reference to the input vertex.
         * @return input vertex
         */
        std::shared_ptr<data_structures::Vertex> getInput();

        /**
         * Get the reference to the output vertex.
         * @return output vertex
         */
        std::shared_ptr<data_structures::Vertex> getOutput();

        /**
         * Get the weight.
         * @return weight
         */
        double getWeight();

        /**
         * Get the number of times this edge can be traversed.
         * @return traverse limit
         */
        unsigned int getTraverseLimit();

        /**
         * Combines the current weight with the new one. By default, the operation is addition.
         * x*w1 + x*w2 = y
         * x * (w1 + w2) = y
         * @param argWeight weight to combine
         */
        void combineWeight(double argWeight);

        /**
         * Gets the value of the input vertex, multiplies is wit the weight, and combines the result with the output
         * vertex.
         * @return result
         */
        double propagateValue();

        /**
         * Has this edge been traversed at least once?
         * @return true, if this edge has been traversed at least once
         */
        bool isTraversed();

        /**
         * Has this edge been traversed the maximum number of allowed times?
         * @return true, if traverseCount equals of is greater than traverseLimit
         */
        bool isAtTaversalLimit();

        /**
         * Sets the traverseCount to 0.
         */
        void reset();

        unsigned int getIndex() const;

        void setIndex(unsigned int index);

        /**
         * Get a string with useful information about this object.
         * @return this object's information
         */
        std::string toString();

        bool isFlaggedForDeletion() const;

        void setFlaggedForDeletion(bool flaggedForDeletion);
    };
}


#endif //NEUROEVOLUTION_EDGE_H
