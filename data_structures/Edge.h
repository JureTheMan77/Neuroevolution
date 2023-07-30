//
// Created by jure on 7/13/20.
//

#ifndef NEUROEVOLUTION_EDGE_H
#define NEUROEVOLUTION_EDGE_H

#include <memory>
#include "Edge.fwd.h"
#include "Vertex.fwd.h"
#include "ICrossoverable.h"
#include "IDeepCloneable.h"

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
         * New edge constructor.
         * @param input input vertex
         * @param output output vertex
         * @param index edge index/unique id
         * @param weight weight of the edge
         * @param traverseLimit how many times this edge can be traversed
         * @param dominant is this edge dominant
         */
        Edge(const std::shared_ptr<data_structures::Vertex> &input,
             const std::shared_ptr<data_structures::Vertex> &output, unsigned int index, double weight,
             unsigned int traverseLimit, bool dominant) : ICrossoverable(dominant), input(input), output(output),
                                                          index(index), weight(weight), traverseLimit(traverseLimit) {};

        /**
         * Creates a new edge shared pointer.
         * @param input input vertex
         * @param output output vertex
         * @param index edge index/unique id
         * @param weight weight of the edge
         * @param traverseLimit how many times this edge can be traversed
         * @param dominant is this edge dominant
         * @return new edge shared pointer
         */
        static std::shared_ptr<data_structures::Edge> createEdge(const std::shared_ptr<data_structures::Vertex> &input,
                                                                 const std::shared_ptr<data_structures::Vertex> &output,
                                                                 unsigned int index, double weight,
                                                                 unsigned int traverseLimit, bool dominant);

        /**
         * Destroys this object.
         */
        ~Edge() = default;

        /**
         * Get the reference to the input vertex. Calls vertex.lock() internally.
         * @return input vertex shared pointer
         */
        [[nodiscard]] std::shared_ptr<data_structures::Vertex> getInput() const;

        /**
         * Get the reference to the output vertex. Calls vertex.lock() internally.
         * @return output vertex shared pointer
         */
        [[nodiscard]] std::shared_ptr<data_structures::Vertex> getOutput() const;

        /**
         * Get the weight.
         * @return weight
         */
        [[nodiscard]] double getWeight() const;

        /**
         * Set the weight.
         * @param weightArg weight so set
         */
        void setWeight(double weightArg);

        /**
         * Get the number of times this edge can be traversed.
         * @return traverse limit
         */
        [[nodiscard]] unsigned int getTraverseLimit() const;

        /**
         * Set the number of times this edge can be traversed.
         * @param traverseLimitArg traverse limit
         */
        void setTraverseLimit(unsigned int traverseLimitArg);

        /**
         * Combines the current weight with the new one. By default, the operation is addition.
         * @param argWeight weight to combine
         */
        void combineWeight(double argWeight);

        /**
         * Gets the value of the input vertex, multiplies is with the weight, and combines the result with the output
         * vertex. If the edge is at it's traverse limit, then don't propagate and return 0.
         * @return result
         */
        double propagateValue();

        /**
         * Has this edge been traversed at least once?
         * @return true, if this edge has been traversed at least once
         */
        [[nodiscard]] bool isTraversedOnce() const;

        /**
         * Has this edge been traversed the maximum number of allowed times?
         * @return true, if traverseCount equals of is greater or equal to traverseLimit
         */
        [[nodiscard]] bool isAtTraverseLimit() const;

        /**
         * How many times can this edge still be traversed?
         * @return traversals remaining
         */
        unsigned int traversalsRemaining() const;

        /**
         * Sets the traverseCount to 0.
         */
        void reset();

        /**
         * Get the index of this edge.
         * @return index
         */
        [[nodiscard]] unsigned int getIndex() const;

        void setIndex(unsigned int indexArg);

        /**
         * Get a string with useful information about this object.
         * @param technical format for https://csacademy.com/app/graph_editor/
         * @return this object's information
         */
        [[nodiscard]] std::string toString(bool technical) const;

        /**
         * Needed for agent minimization step.
         * @return should this edge be kept
         */
        [[nodiscard]] bool isFlaggedForDeletion() const;

        /**
         * Needed for agent minimization step.
         * @param flaggedForDeletionArg should this edge be kept
         */
        void setFlaggedForDeletion(bool flaggedForDeletionArg);
    };
}


#endif //NEUROEVOLUTION_EDGE_H
